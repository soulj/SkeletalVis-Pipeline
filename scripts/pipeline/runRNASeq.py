#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:25:47 2017

@author: Jamie Soul

Script used to drive Galaxy to analyse RNA-Seq data.

Need to define a csv with headers accession number, species, foldChangeOnly, stranded and single.
Analysis will be done once for each accession number. Other args are used as input to tools in
the pipeline. Also need to define a comparisons file and sample tables. These are used as
input to the tools and to identify what sample particular files belong to.

We first get the data from the Galaxy repository and pass it through fastqc tool for a quality read.
We then run it through suite of kallisto tools to extract a gene level counts matrix. We also 
calculate the fold change and then cut the produced table to get a fold change for each
comparison.

The cut tables then get passed one by one into the downstream pathway analysis workflow

To run you need to set an api_key, see https://galaxyproject.org/develop/api/
for details. You also need to define a csv of inputs, a comparisons file and a sample table.
"""
def runWorkflow(argDictionary, comparisons,samples):
    from bioblend.galaxy import GalaxyInstance
    from bioblend.galaxy.histories import HistoryClient
    from bioblend.galaxy.tools import ToolClient
    from bioblend.galaxy.workflows import WorkflowClient
    from bioblend.galaxy.libraries import LibraryClient
    import tempfile
    
    
    import time
    api_key = ''
    galaxy_host = ''

    gi = GalaxyInstance(url=galaxy_host, key=api_key)

    history_client = HistoryClient(gi)
    tool_client = ToolClient(gi)
    workflow_client = WorkflowClient(gi)
    library_client = LibraryClient(gi)
    
    history = history_client.create_history(argDictionary['accessionNumber'])
    
    comparisonsTable = tool_client.upload_file(comparisons, history['id'], file_type='txt')
    sampleTable = tool_client.upload_file(samples, history['id'], file_type='tabular')
    
    if argDictionary['site'] == "ENA":
        #fastqs available on ENA    
        tool_inputs = {
                "accessionNumber":argDictionary["ENA"],"sampleTable":{'id': sampleTable['outputs'][0]['id'], 'src': 'hda'}
                
            }
        
    
        #run the tool to get the data from ENA
        tool_client.run_tool(history['id'],'getRNASeqExpressionData', tool_inputs)
        
        #we want to wait until we have all datasets
        while getNumberNotComplete(history['id'], history_client) > 0:
            time.sleep(10)
            
        
        #sleep until all the fastq files are findable
        time.sleep(120)
        
        
        dirpath = tempfile.mkdtemp()
        fileList = getDatasetsByApproxName("files.tabular", history,history_client)[0]
        fileList = history_client.download_dataset(history["id"],fileList["id"],dirpath)
        num_lines = sum(1 for line in open(fileList)) -1
        
        datasets=list()
        while len(datasets)!=num_lines:
                    time.sleep(10)
                    datasets = getDatasetsByApproxName("fastq",history,history_client )                
    else: #for SRA       
    
        if argDictionary['single'] == "TRUE":
            with open(samples) as tsvfile:
                reader = csv.DictReader(tsvfile, delimiter='\t')
                for sample in reader:
                    print (sample)
                    fileNames=str.split(sample["File"],"|")
                    for fileName in fileNames:                    
                        tool_inputs = {
                                "input|input_select":"accession_number",
                                "outputformat":"fastqsanger.gz",
                                "input|accession":fileName   
                            }
                        #run the tool to get the single data from SRA
                        tool_client.run_tool(history['id'],'toolshed.g2.bx.psu.edu/repos/iuc/sra_tools/fastq_dump/2.8.1.3', tool_inputs)
               
        else:
             with open(samples) as tsvfile:
                reader = csv.DictReader(tsvfile, delimiter='\t')
           
                for sample in reader:            
                    tool_inputs = {
                            "accession_number":sample["File"]           
                        }
                    #run the tool to get the paired data from SRA
                    tool_client.run_tool(history['id'],'toolshed.g2.bx.psu.edu/repos/mandorodriguez/fastqdump_paired/fastq_dump_paired/1.1.4', tool_inputs)
                
        while getNumberNotComplete(history['id'], history_client) > 0:
            time.sleep(10)
     
    datasets = getDatasetsByApproxName("fastq",history,history_client )
    #get the fastQC tool
    for fastq in datasets:
        try:
            tool_inputs = {'input_file' : {'id': fastq['id'], 'src': 'hda'}}
            tool_client.run_tool(history['id'],'toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.69', tool_inputs)
        except Exception:
            pass
        
    #wait till complete
    while getNumberNotComplete(history['id'], history_client) > 0:
        time.sleep(10)
    
    #make dataset collections for quantification using the fastq files
    collections=list()
    with open(samples) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            datasets=list()
            fileNames=str.split(row["File"],"|")
            
            for fileName in fileNames:
                datasets= datasets + getDatasetsByApproxName(fileName,history,history_client )
                    
            #make list of datasets
            collections.append(makeDataSetCollection(datasets,row["Sample"],history,history_client))
            
            
            
    #get the correct kallisto index
    species = argDictionary['species'].lower()
    index = getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name=species +"IndexFile")
    index = {'id': index, 'src': 'hda'}
    
    #run kallisto for every dataset collection
    for collection in collections:
        #set up the tool_inputs
        tool_inputs = {'index' : index,'inputs' : {'id': collection['id'], 'src': 'hdca'} ,"single":argDictionary["single"],"stranded":argDictionary["stranded"]}
        
        
        #often encounter connection broken error - possible problem with Certus server?
        #bypass by ignoring the exception
        tool_client.run_tool(history['id'],'kallistoQuant', tool_inputs)


    # we want to wait until we have all datasets
    while getNumberNotComplete(history['id'], history_client) > 0:
        time.sleep(10)
        
    # Run multiqc on kallisto logs and fastqc files
    datasets = getDatasetsByApproxName("RawData",history,history_client )
    kallistoLogs = getDatasetsByApproxName(".log", history, history_client)
    
    tool_inputs = {}
    for i, dataset in enumerate(datasets+kallistoLogs):
        if not dataset["deleted"]:
            if dataset in datasets:
                software = 'fastqc'
            else:
                software = 'kallisto'
            params = {'id' : dataset['id'], 'src': 'hda', 'name': dataset['name']}
            tool_inputs.update({'results_%s|software_cond|software' % i: software, 'results_%s|input_file' % i: params})

#    #summarise with the multiQC tool
    tool_client.run_tool(history['id'],'multiqc', tool_inputs)
    
    multiQc = getDatasetsByApproxName("multiqc",history,history_client)[0]
    
        
    #get all the abundance files to convert to gene level counts matrix
    datasets = getDatasetsByApproxName(".abundance",history,history_client )
    
    #make a dataset collection for to make a countsMatrix
    collection = makeDataSetCollection(datasets,"abundances",history,history_client)
    
    
    #set up the tool_inputs
    tool_inputs = {'inputs' : {'id': collection['id'], 'src': 'hdca'} ,"species":argDictionary['species']}
    
    #convert abundances to gene level counts matrix
    tool_client.run_tool(history['id'],'KallistoAbundancestoGeneCountMatrix', tool_inputs)
    
    # A diry hack, we want to wait until we have all datasets
    while getNumberNotComplete(history['id'], history_client) > 0:
        time.sleep(10)
    
    txi = getDatasetsByApproxName("txi",history,history_client)
    

    #set up the tool_inputs for PCA
    tool_inputs = {'txiData' : {'id': txi[0]['id'], 'src': 'hda'} ,'sampleTable' : {'id': sampleTable['outputs'][0]['id'], 'src': 'hda'} ,"species":argDictionary['species'],'technicalReplicates':argDictionary['technicalReplicates'],'batchCorrect':argDictionary['batchCorrect']}
    
    #run deseq2
    tool_client.run_tool(history['id'],'PCARNASeq', tool_inputs)
    
    pca = getDatasetsByApproxName("PCA",history,history_client)[0]
    
       
    #set up the tool_inputs for DESeq2
    tool_inputs = {'txiData' : {'id': txi[0]['id'], 'src': 'hda'} ,'sampleTable' : {'id': sampleTable['outputs'][0]['id'], 'src': 'hda'} ,
    'comparisonsTable' : {'id': comparisonsTable['outputs'][0]['id'], 'src': 'hda'} ,"foldChangeOnly":argDictionary['foldChangeOnly'],"species":argDictionary['species'],'technicalReplicates':argDictionary['technicalReplicates'],'batchCorrect':argDictionary['batchCorrect']}
    
    #run deseq2
    tool_client.run_tool(history['id'],'DESeq2FoldChange', tool_inputs)
         
    #run chrdir
    tool_client.run_tool(history['id'],'characteristicDirectionRNASeq', tool_inputs)
    
        #we want to wait until we have all datasets
    while getNumberNotComplete(history['id'], history_client) > 0:
        time.sleep(10)
        
        
    #get the foldchange data, cut and run pathway workflow    
    dataset_id = getFoldChangeData(history, history_client)['id']
    
    
    return_collection = [{'accessionNo':argDictionary['accessionNumber'], 'foldChange': getUrl(dataset_id), 'PCA': getUrl(pca["id"]),'chrDirTable': getUrl(getMostRecentDatasetByName('chrDirTable.tabular', history, history_client)['id'])}]
    
    
    number_of_comparisons = -1
    for line in open(comparisons):
        if not line.isspace():
            number_of_comparisons += 1

    for comparison in range(0, int(number_of_comparisons)):
        tool_inputs = {
            'foldChangeTable' : {'id': dataset_id, 'src': 'hda'},
            'comparisonNumber' : comparison + 1
        }
        tool_client.run_tool(history['id'], 'cutFoldChangeTable', tool_inputs)
        
    while getNumberNotComplete(history['id'], history_client) > 0:
        time.sleep(10)
        
        
    if argDictionary['species'] in ["Rat","Cow","Horse","Pig","Zebrafish"]:
        pathwayAnalysisWorkflow = workflow_client.show_workflow('c9468fdb6dc5c5f1')
        
        params = dict()
        for key in pathwayAnalysisWorkflow['steps'].keys():
            params[key] = argDictionary
        
        if argDictionary['species'] == "Rat":
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="ratStringNetwork")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="ratGeneLengths")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="HOM_AllOrganism.rpt")
        if argDictionary['species'] == "Cow":
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="cowStringNetwork")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="cowGeneLengths")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="HOM_AllOrganism.rpt")
        if argDictionary['species'] == "Horse":
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="horseStringNetwork")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="horseGeneLengths")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="Homology.horse.txt")
        if argDictionary['species'] == "Pig":
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="pigStringNetwork.txt")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="pigGeneLengths.tabular")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="Homology.pig.txt")
        if argDictionary['species'] == "Zebrafish":
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="zebrafishStringNetwork")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="zebrafishGeneLengths")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="HOM_AllOrganism.rpt")
        
                
        pathwayDatamap = {'3' : {'id': homology, 'src': 'hda'},'2' : {'id': network, 'src': 'hda'},'1' : {'id': geneLengths, 'src': 'hda'}}

        diffExpDataCollection = getDatasetsByName('cutTable.tabular', history, history_client)
        for index, diffExpData in enumerate(diffExpDataCollection):
            
            numCompleted = getNumberComplete(history['id'], history_client) + 10
            print(numCompleted)
            
            pathwayDatamap["0"] = {'id': diffExpData['id'], 'src': 'hda'}
            workflow_client.invoke_workflow(pathwayAnalysisWorkflow['id'], 
                                            inputs = pathwayDatamap, 
                                            history_id = history['id'], 
                                            params = params)                  
            
            
            comparisonDict = getRowFromCsv(comparisons, index)
            
            if 'Factor1' in comparisonDict.keys():
                comparisonDict['Factor'] = comparisonDict['Factor1'] + "." + comparisonDict['Factor2']
                
            return_dict = {'accessionNo':argDictionary['accessionNumber'],
                           'factor':comparisonDict['Factor'],
                           'comparisonNum':comparisonDict['Numerator'],
                           'comparisonDenom':comparisonDict['Denominator'],
                           'foldChange': getUrl(diffExpData['id']),
                           'interactome': pathwayDatamap['0']['id'],
                           'exonLength': pathwayDatamap['2']['id']}
            
            while getNumberComplete(history['id'], history_client) < numCompleted:
                time.sleep(10)
    
            return_dict['moduleNodes'] = getUrl(getMostRecentDatasetByName('moduleNodes.text', 
                history, history_client)['id'])
            return_dict['modulePlots'] = getUrl(getMostRecentDatasetByName('modulePlots.pdf',
                history, history_client)['id'])
            return_dict['slimEnrichPathways'] = getUrl(getMostRecentDatasetByName('slimEnrichmentPathways.tabular',
                history, history_client)['id'])
            return_dict['enrichedDrugsReverse'] = getUrl(getMostRecentDatasetByName('enrichedDrugsReverse.tabular',
                history, history_client)['id'])
            return_dict['enrichedDrugsMimic'] = getUrl(getMostRecentDatasetByName('enrichedDrugsMimic.tabular',
                history, history_client)['id'])
            return_dict['enrichedTerms'] = getUrl(getMostRecentDatasetByName('enrichedTerms.tabular',
                history, history_client)['id'])
            return_dict['enrichedTerms.reduced'] = getUrl(getMostRecentDatasetByName('enrichedTerms.reduced.tabular',
                history, history_client)['id'])
            return_dict['GO.MDS'] = getUrl(getMostRecentDatasetByName('GO.MDS.html',
                history, history_client)['id'])
            return_collection.append(return_dict)
       
        # Hard code keys to define the order
        keys = ['accessionNo','multiQC','factor','PCA','chrDirTable','comparisonNum','comparisonDenom','foldChange',
        'interactome','exonLength','moduleNodes','modulePlots',
        'slimEnrichPathways','secretedProteins','enrichedDrugsReverse','enrichedDrugsMimic','enrichedTerms','enrichedTerms.reduced','GO.MDS']
        
        outFileName = 'output/' +  argDictionary['accessionNumber'] + '-workflowOutput.tsv'
        
        with open(outFileName, 'wb') as csvFile:
            # Get headers from last dictionary in collection as first doesn't contain all keys
            csvOutput = csv.DictWriter(csvFile, keys, delimiter = "\t")
            csvOutput.writeheader()
            csvOutput.writerows(return_collection)
            
        #tool_client.upload_file(outFileName, history['id'], file_type='tsv')
        
        return return_collection
    else:  
        pathwayAnalysisWorkflow = workflow_client.show_workflow('e85a3be143d5905b')
        
        params = dict()
        for key in pathwayAnalysisWorkflow['steps'].keys():
            params[key] = argDictionary
            
       
        if argDictionary['species'] == "Mouse":  
        
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="mouseStringNetwork")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="MouseGeneLengths.tab")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="Homology.mouse.txt")
            secretedReference=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="uniprot-secreted-mouse.txt")
            
            pathwayDatamap = {'4' : {'id':  secretedReference, 'src': 'hda'},'3' : {'id': homology, 'src': 'hda'},'2' : {'id': network, 'src': 'hda'},'1' : {'id': geneLengths, 'src': 'hda'}}
        else:
        
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="humanStringNetwork")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="geneLengths")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="Homology.mouse.txt")
            secretedReference=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="uniprot-secreted-human.txt")
            pathwayDatamap = {'4' : {'id':  secretedReference, 'src': 'hda'},'3' : {'id': homology, 'src': 'hda'},'2' : {'id': network, 'src': 'hda'},'1' : {'id': geneLengths, 'src': 'hda'}}
    
        diffExpDataCollection = getDatasetsByName('cutTable.tabular', history, history_client)
        for index, diffExpData in enumerate(diffExpDataCollection):
            
            numCompleted = getNumberComplete(history['id'], history_client) + 14
            print(numCompleted)
            
            pathwayDatamap["0"] = {'id': diffExpData['id'], 'src': 'hda'}

    
        
            #pathwayDatamap['1'] = {'id': diffExpData['id'], 'src': 'hda'}
            workflow_client.invoke_workflow(pathwayAnalysisWorkflow['id'], 
                                            inputs = pathwayDatamap, 
                                            history_id = history['id'], 
                                            params = params)
            comparisonDict = getRowFromCsv(comparisons, index)
            
            if 'Factor1' in comparisonDict.keys():
                comparisonDict['Factor'] = comparisonDict['Factor1'] + "." + comparisonDict['Factor2']
                
            return_dict = {'accessionNo':argDictionary['accessionNumber'],
                           'factor':comparisonDict['Factor'],
                           'comparisonNum':comparisonDict['Numerator'],
                           'comparisonDenom':comparisonDict['Denominator'],
                           'foldChange': getUrl(diffExpData['id']),
                           'interactome': pathwayDatamap['0']['id'],
                           'exonLength': pathwayDatamap['2']['id']}
            
            while getNumberComplete(history['id'], history_client) < numCompleted:
                time.sleep(10)
    
            return_dict['moduleNodes'] = getUrl(getMostRecentDatasetByName('moduleNodes.text', 
                history, history_client)['id'])
            return_dict['modulePlots'] = getUrl(getMostRecentDatasetByName('modulePlots.pdf',
                history, history_client)['id'])
            return_dict['pathways'] = getUrl(getMostRecentDatasetByName('pathways.tabular', 
                history, history_client)['id'])
            return_dict['enrichPlot'] = getUrl(getMostRecentDatasetByName('enrichmentPlot.png', 
                history, history_client)['id'])
            return_dict['enrichmentTable'] = getUrl(getMostRecentDatasetByName('TF_EnrichmentTable.tabular', 
                history, history_client)['id'])
            return_dict['slimEnrichPathways'] = getUrl(getMostRecentDatasetByName('slimEnrichmentPathways.tabular',
                history, history_client)['id'])
            return_dict['secretedProteins'] = getUrl(getMostRecentDatasetByName('secretedProteins.tabular',
                history, history_client)['id'])
            return_dict['enrichedDrugsReverse'] = getUrl(getMostRecentDatasetByName('enrichedDrugsReverse.tabular',
                history, history_client)['id'])
            return_dict['enrichedDrugsMimic'] = getUrl(getMostRecentDatasetByName('enrichedDrugsMimic.tabular',
                history, history_client)['id'])
            return_dict['enrichedTerms'] = getUrl(getMostRecentDatasetByName('enrichedTerms.tabular',
                history, history_client)['id'])
            return_dict['enrichedTerms.reduced'] = getUrl(getMostRecentDatasetByName('enrichedTerms.reduced.tabular',
                history, history_client)['id'])
            return_dict['GO.MDS'] = getUrl(getMostRecentDatasetByName('GO.MDS.html',
                history, history_client)['id'])
            return_collection.append(return_dict)
       
        # Hard code keys to define the order
        keys = ['accessionNo','multiQC','factor','PCA','chrDirTable','comparisonNum','comparisonDenom','foldChange',
        'interactome','exonLength','moduleNodes','modulePlots','pathways','enrichPlot', 'enrichmentTable',
        'slimEnrichPathways','secretedProteins','enrichedDrugsReverse','enrichedDrugsMimic','enrichedTerms','enrichedTerms.reduced','GO.MDS']
        
        outFileName = 'output/' +  argDictionary['accessionNumber'] + '-workflowOutput.tsv'
        
        with open(outFileName, 'wb') as csvFile:
            # Get headers from last dictionary in collection as first doesn't contain all keys
            csvOutput = csv.DictWriter(csvFile, keys, delimiter = "\t")
            csvOutput.writeheader()
            csvOutput.writerows(return_collection)
            
        
        return return_collection

def getFoldChangeData(history, history_client):
    return getMostRecentDatasetByName('foldChangeTable.tabular', history, history_client)


def getLibraryData(history,samples,history_client,library_client):
    
    datasets = library_client.show_library(contents=True,library_id='2d9035b3fc152403')
    
    with open(samples) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            for dataset in datasets:
                if row["File"] in dataset["name"]:
                   if "unmapped" not in dataset["name"]:
                       if ".bam" not in dataset["name"]:
                           if "alignment" not in dataset["name"]:
                               history_client.upload_dataset_from_library(history["id"],dataset["id"])
    
def getLibraryToolDataID(library_client,history_client,history,name):
    
    datasets = library_client.show_library(contents=True,library_id='c9468fdb6dc5c5f1')
    
    for dataset in datasets:
        if name in dataset["name"]:
           history_client.upload_dataset_from_library(history["id"],dataset["id"])
           dataset=getMostRecentDatasetByName(dataset["name"].split("/")[-1],history,history_client)
           return dataset["id"]



def getDatasetsByName(datasetName, history, history_client):
    datasets = history_client.show_history(history['id'], contents = True)
    matchingDatasets = []
    for dataset in datasets:
        if datasetName == dataset['name']:
             if not dataset["deleted"]:
                 matchingDatasets.append(dataset)
    return matchingDatasets   

def getDatasetsByApproxName(datasetName, history, history_client):
    datasets = history_client.show_history(history['id'], contents = True)
    matchingDatasets = []
    for dataset in datasets:
        if datasetName in dataset['name']:
            if not dataset["deleted"]:
                matchingDatasets.append(dataset)
    return matchingDatasets  

   
def getMostRecentDatasetByName(datasetName, history, history_client):
    matchingDatasets = getDatasetsByName(datasetName, history, history_client)
    if len(matchingDatasets) == 1:
        return matchingDatasets[0]
    elif len(matchingDatasets) == 0:
        return {}
    else:
        return matchingDatasets[-1]    

def getMostRecentDataset(dataCollection,history_client):
    mostRecentData = dataCollection[0]
    mostRecentDataFull = history_client.show_dataset(mostRecentData['id'],mostRecentData['history_id'])
    for data in dataCollection:
        dataFull = history_client.show_dataset(data['id'],data['history_id'])
        if dataFull['create_time'] > mostRecentDataFull['create_time']:
            mostRecentData = data
    return mostRecentData


def makeDataSetCollection(datasets,name,history,history_client):
    
    dataList = list()
    for dataset in datasets:
        if not dataset["deleted"]:
            dataList.append({'id': dataset["id"],'name': dataset["name"],'src': 'hda'})

    collection = {'collection_type':'list','element_identifiers' :  dataList,'name' : name }    
    collection = history_client.create_dataset_collection(history['id'],collection)
    
   
    return collection
    
        

# Need a way to delay before continuing otherwise results in an error due to missing 'id' as the
# dataset hasn't been created yet. This could be a bit of an overly strong solution, the
# processing will only continue when everything has completed. 
def getNumberNotComplete(historyId, history_client):
    try:
        states = history_client.show_history(historyId)['state_ids']
        return len(states['running']) + len(states['queued']) + len(states['new']) + len(states['paused'])
    except Exception:
        return 1
    

def getNumberComplete(historyId, history_client):
    states = history_client.show_history(historyId)['state_ids']
    return len(states['ok'])

def getUrl(identifier):
    return 'http://localhost:8080/' + identifier + '/display/?preview=True'
    
   
def getRowFromCsv(csvPath, rowNo):
    import itertools
    
    with open(csvPath) as csvfile:
        for rw in itertools.islice(csv.DictReader(csvfile, delimiter='\t'), rowNo, rowNo+1):
            return rw
    

import csv

with open('exampleRNASeqParams.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        hst = runWorkflow(row, 'comparisons/' + row['accessionNumber'] + '_comparisons.txt','sampleTables/' + row['accessionNumber'] + '_sampleTable.txt')
        


