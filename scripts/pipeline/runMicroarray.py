# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 16:57:17 2017

@author: Robert Ashton

Script used to run two workflows, the first 'Differential Expression Workflow' 
downloads data from ArrayExpress or GEO and run differential expression tool to
get table of fold change and p values.

Then cut table into separate comparisons and run 'Downstream pathway analysis'
workflow on each one. Then writes URLs of output to a csv file.

To run you need to set an api_key, see https://galaxyproject.org/develop/api/
for details. You also need to import copies of the workflows, these can be
found in the 'Shared Data' section of Galaxy. You then need to use the ID of
these as arguments to workflow_client.show_workflow.

"""
def runWorkflow(argDictionary, comparisons):
    from bioblend.galaxy import GalaxyInstance
    from bioblend.galaxy.histories import HistoryClient
    from bioblend.galaxy.tools import ToolClient
    from bioblend.galaxy.workflows import WorkflowClient
    from bioblend.galaxy.libraries import LibraryClient
    import time
    
    api_key = ''
    galaxy_host = 'http://localhost:8080/'

    gi = GalaxyInstance(url=galaxy_host, key=api_key)

    history_client = HistoryClient(gi)
    tool_client = ToolClient(gi)
    workflow_client = WorkflowClient(gi)
    library_client = LibraryClient(gi)
    
    history = history_client.create_history(row['accessionNumber'])
    # Import the galaxy workflow
    workflow = workflow_client.show_workflow('a799d38679e985db')

    input_file = tool_client.upload_file(comparisons, history['id'], file_type='txt')

    # Run workflow on csv data to create a new history.
    params = dict()
    for key in workflow['steps'].keys():
        params[key] = argDictionary
    
    datamap = {'1' : {'id': input_file['outputs'][0]['id'], 'src': 'hda'}}

    workflow_client.invoke_workflow(workflow['id'], inputs = datamap, history_id = history['id'], params = params)
    
    # A diry hack, we want to wait until we have all datasets
    while getNumberNotComplete(history['id'], history_client) > 0:
        time.sleep(10)
        
    
    dataset_id = getFoldChangeData(history, history_client)['id']

    
    return_collection = [{'accessionNo':argDictionary['accessionNumber'], 'foldChange': getUrl(dataset_id),
    'PCA': getUrl(getMostRecentDatasetByName('PCAplot.png', history, history_client)['id']),'chrDirTable': getUrl(getMostRecentDatasetByName('chrDirTable.tabular', history, history_client)['id'])}]
    
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
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="Homology.rat.txt")
        if argDictionary['species'] == "Cow":
            network=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="cowStringNetwork")
            geneLengths=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="cowGeneLengths")
            homology=getLibraryToolDataID(history=history,history_client=history_client,library_client=library_client,name="Homology.cow.txt")
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
                
            if 'Paired1' in comparisonDict.keys():
                comparisonDict['Factor'] = comparisonDict['Paired1']
                
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
            return_dict['slimEnrichmentPathways'] = getUrl(getMostRecentDatasetByName('slimEnrichmentPathways.tabular',
            history, history_client)['id'])
            return_dict['slimEnrichmentPlot'] = getUrl(getMostRecentDatasetByName('slimEnrichmentPlot.png',
            history, history_client)['id'])
            return_collection.append(return_dict)     
       
        # Hard code keys to define the order
        keys = ['accessionNo','factor','comparisonNum','comparisonDenom','PCA','chrDirTable','foldChange',
        'interactome','exonLength','moduleNodes','modulePlots','enrichmentTable','slimEnrichmentPathways','slimEnrichmentPlot']
        with open('output/' +  argDictionary['accessionNumber'] + '-workflowOutput.csv', 'wb') as csvFile:
            # Get headers from last dictionary in collection as first doesn't contain all keys
            csvOutput = csv.DictWriter(csvFile, keys)
            csvOutput.writeheader()
            csvOutput.writerows(return_collection)
            
        return return_collection
    else: 
        pathwayAnalysisWorkflow = workflow_client.show_workflow('e85a3be143d5905b')
        
        params = dict()
        for key in pathwayAnalysisWorkflow['steps'].keys():
            params[key] = argDictionary
            
        # MouseGeneLengths.tab has id 457f69dd7016f307 - step 2 of workflow
        # Mouse interactome has id 073be90ac6c3bce5 - step 0 of workflow
        
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

            workflow_client.invoke_workflow(pathwayAnalysisWorkflow['id'], 
                                            inputs = pathwayDatamap, 
                                            history_id = history['id'], 
                                            params = params)                  
            
            
            comparisonDict = getRowFromCsv(comparisons, index)
            
            if 'Factor1' in comparisonDict.keys():
                comparisonDict['Factor'] = comparisonDict['Factor1'] + "." + comparisonDict['Factor2']
                
            if 'Paired1' in comparisonDict.keys():
                comparisonDict['Factor'] = comparisonDict['Paired1']
                
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
            return_dict['slimEnrichmentPathways'] = getUrl(getMostRecentDatasetByName('slimEnrichmentPathways.tabular',
            history, history_client)['id'])
            return_dict['slimEnrichmentPlot'] = getUrl(getMostRecentDatasetByName('slimEnrichmentPlot.png',
            history, history_client)['id'])
            return_collection.append(return_dict)     
       
        # Hard code keys to define the order
        keys = ['accessionNo','factor','comparisonNum','comparisonDenom','PCA','chrDirTable','foldChange',
        'interactome','exonLength','moduleNodes','modulePlots','pathways','enrichPlot','enrichmentTable','slimEnrichmentPathways','slimEnrichmentPlot']
        with open('output/' +  argDictionary['accessionNumber'] + '-workflowOutput.csv', 'wb') as csvFile:
            # Get headers from last dictionary in collection as first doesn't contain all keys
            csvOutput = csv.DictWriter(csvFile, keys)
            csvOutput.writeheader()
            csvOutput.writerows(return_collection)
            
        return return_collection

def getFoldChangeData(history, history_client):
    return getMostRecentDatasetByName('foldChangeTable.tabular', history, history_client)


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
    
# Need a way to delay before continuing otherwise results in an error due to missing 'id' as the
# dataset hasn't been created yet. This could be a bit of an overly strong solution, the
# processing will only continue when everything has completed. 
def getNumberNotComplete(historyId, history_client):
    states = history_client.show_history(historyId)['state_ids']
    return len(states['running']) + len(states['queued']) + len(states['new']) + len(states['paused'])

def getNumberComplete(historyId, history_client):
    states = history_client.show_history(historyId)['state_ids']
    return len(states['ok'])

def getUrl(identifier):
    return 'http://localhost:8080/datasets/' + identifier + '/display/?preview=True'
    
def getRowFromCsv(csvPath, rowNo):
    import itertools
    
    with open(csvPath) as csvfile:
        for rw in itertools.islice(csv.DictReader(csvfile, delimiter='\t'), rowNo, rowNo+1):
            return rw
    

import csv

with open('exampleMicroarrayParams.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        hst = runWorkflow(row, 'comparisons/' + row['accessionNumber'] + '_comparisons.txt')
