#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 13:24:09 2017

Given an accession number - select the most recent history - get the experiment level data like QC and PCA
then for the comparisons get the fold changes, pathway, networks etc. Number _1 _2 etc


@author: Jamie Soul
"""

from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.histories import HistoryClient
import urllib
import os
import zipfile
import csv
import itertools


def getDatasetsByName(datasetName, history, history_client):
    datasets = history_client.show_history(history['id'], contents = True)
    matchingDatasets = []
    for dataset in datasets:
        if datasetName == dataset['name']:
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


def getPCA(history,history_client,directory,galaxy_host):
    pca = getDatasetsByApproxName('PCA', history, history_client)[0]
    pcaFile = urllib.URLopener()
    url = galaxy_host + pca["url"] + "/display?to_ext=png"
    out = directory + "/" + history["name"] + "_PCA.png"
    pcaFile.retrieve(url,out)

def getChrDirTable(history,history_client,directory,galaxy_host):
    pca = getDatasetsByApproxName('chrDirTable.tabular', history, history_client)[0]
    pcaFile = urllib.URLopener()
    url = galaxy_host + pca["url"] + "/display?to_ext=html"
    out = directory + "/" + history["name"] + "_chrDirTable.txt"
    pcaFile.retrieve(url,out)

def getQC(history,history_client,directory,galaxy_host):
    qc = getDatasetsByName('microarrayQC.html', history, history_client)[0]
    qcFile = urllib.URLopener()
    url = galaxy_host + qc["url"] + "/display?to_ext=html"
    out = directory + "/" + history["name"] + "_microarrayQC.zip"
    qcFile.retrieve(url,out)

    zip_ref = zipfile.ZipFile(out, 'r')
    zip_ref.extractall(directory)
    zip_ref.close()
    
    
def getComparisonsTable(history,history_client,directory,galaxy_host):
    comp = getDatasetsByApproxName('comparisons', history, history_client)[0]
    compFile = urllib.URLopener()
    url = galaxy_host + comp["url"] + "/display?to_ext=html"
    out = directory + "/" + history["name"] + "_comparisons.txt"
    compFile.retrieve(url,out)
    return(out)

def getFoldChange(i,history,history_client,dataDirectory,galaxy_host):
    fc = getDatasetsByApproxName('cutTable', history, history_client)[i]
    fcFile = urllib.URLopener()
    url = galaxy_host + fc["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_FC_" + str(i+1) + ".txt"
    fcFile.retrieve(url,out)

def getPathways(i,comparison,pvalue,foldchange,history,history_client,dataDirectory,galaxy_host):
    pathways = getDatasetsByApproxName('slimEnrichmentPathways', history, history_client)[i]
    pathwaysFile = urllib.URLopener()
    url = galaxy_host +  pathways["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_pathways_" + str(comparison+1) + "_" + foldchange + "_" + pvalue + ".txt"
    pathwaysFile.retrieve(url,out)
    
def getTFs(i,comparison,pvalue,foldchange,history,history_client,dataDirectory,galaxy_host):
    TFs = getDatasetsByApproxName('TF_EnrichmentTable', history, history_client)[i]
    TFsFile = urllib.URLopener()
    url = galaxy_host +  TFs["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_TFs_" + str(comparison+1) + "_" + foldchange + "_" + pvalue +".txt"
    TFsFile.retrieve(url,out)
    
def getStringNetworks(i,comparison,history,history_client,dataDirectory,galaxy_host):
    networks = getDatasetsByName('visNetworks.rdata', history, history_client)[i]
    networksFile = urllib.URLopener()
    url = galaxy_host + networks["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_networks_" + str(comparison+1) +".rdata"
    networksFile.retrieve(url,out)
    networks = getDatasetsByName('summaryTable.text', history, history_client)[i]
    networksFile = urllib.URLopener()
    url = galaxy_host + networks["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_summaryTable_" + str(comparison+1)  +".txt"
    networksFile.retrieve(url,out)

def getBioGridNetworks(i,comparison,history,history_client,dataDirectory,galaxy_host):
    networks = getDatasetsByApproxName('biogridvisNetworks', history, history_client)[i]
    networksFile = urllib.URLopener()
    url = galaxy_host + networks["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_biogridnetworks_" + str(comparison+1) +".rdata"
    networksFile.retrieve(url,out)
    networks = getDatasetsByApproxName('biogridsummaryTable', history, history_client)[i]
    networksFile = urllib.URLopener()
    url = galaxy_host + networks["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_biogridsummaryTable_" + str(comparison+1) +".txt"
    networksFile.retrieve(url,out)

def getDrugEnrichment(i,comparison,pvalue,foldchange,history,history_client,dataDirectory,galaxy_host):
    drugMimic = getDatasetsByApproxName('enrichedDrugsMimic', history, history_client)[i]
    drugReverse = getDatasetsByApproxName('enrichedDrugsReverse', history, history_client)[i]
    drugMimicFile = urllib.URLopener()
    drugReverseFile = urllib.URLopener()
    
    url = galaxy_host + drugMimic["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_drugsMimic_" + str(comparison+1) + "_" + foldchange + "_" + pvalue +".txt"
    drugMimicFile.retrieve(url,out)
    
    url = galaxy_host + drugReverse["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_drugsReverse_" + str(comparison+1) + "_" + foldchange + "_" + pvalue +".txt"
    drugReverseFile.retrieve(url,out)
    
def getGOEnrichment(i,comparison,pvalue,foldchange,history,history_client,dataDirectory,wwwDirectoryPlots,galaxy_host):
    terms = getDatasetsByApproxName('enrichedTerms.tabular', history, history_client)[i]
    termsReduced = getDatasetsByApproxName('enrichedTerms.reduced', history, history_client)[i]
    GOMDS = getDatasetsByApproxName('GO.MDS', history, history_client)[i]
    
    termsFile = urllib.URLopener()
    url = galaxy_host +  terms["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_goterms_" + str(comparison+1) + "_" + foldchange + "_" + pvalue +".txt"
    termsFile.retrieve(url,out)
    
    termsReducedFile = urllib.URLopener()
    url = galaxy_host +  termsReduced["url"] + "/display?to_ext=html"
    out = dataDirectory + "/" + history["name"] + "_goterms_reduced_" + str(comparison+1) + "_" + foldchange + "_" + pvalue + ".txt"
    termsReducedFile.retrieve(url,out)#
    
    GOMDSFile = urllib.URLopener()
    url = galaxy_host +  GOMDS["url"] + "/display?to_ext=html"
    out = wwwDirectoryPlots  + "/" + history["name"] + "_GOMDS_" + str(comparison+1) + "_" + foldchange + "_" + pvalue +".zip"
    GOMDSFile.retrieve(url,out)
    
    wwwDirectoryPlots = wwwDirectoryPlots + "/plot" + str(comparison+1) + "_" + foldchange + "_" + pvalue +"_html"
    zip_ref = zipfile.ZipFile(out, 'r')
    zip_ref.extractall(wwwDirectoryPlots)
    zip_ref.close()


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
    
def getGalaxyData(accession,dataType,species,foldChangeOnly):    
        
    api_key = 'ENTER_API_KEY'
    galaxy_host = 'http://localhost:8080'
    
    gi = GalaxyInstance(url=galaxy_host, key=api_key)
    
    history_client = HistoryClient(gi)

    dataDirectory = "Sybil/Shiny/data/" + accession
    if not os.path.exists(dataDirectory):
        os.makedirs(dataDirectory)
        
    wwwDirectory = "Shiny/www/microarrayQC.html/" + accession
    if not os.path.exists(wwwDirectory):
        os.makedirs(wwwDirectory)
    
    wwwDirectoryPlots = "Shiny/www/plots/" + accession
    if not os.path.exists(wwwDirectoryPlots):
        os.makedirs(wwwDirectoryPlots)
    
    #get the most recent history
    history = history_client.get_histories(name=accession)[0]
    
    #get the experiment level data
    getPCA(history,history_client,dataDirectory,galaxy_host)
    getChrDirTable(history,history_client,dataDirectory,galaxy_host)
    if dataType == "Microarray":
        getQC(history,history_client,wwwDirectory,galaxy_host)
    comparisons = getComparisonsTable(history,history_client,dataDirectory,galaxy_host)
    
    number_of_comparisons = -1
    for line in open(comparisons):
        if not line.isspace():
            number_of_comparisons += 1
            
    if foldChangeOnly=="FALSE":
        pvalues = ["1","0.05"]
        foldchanges = ["1", "1.5","2"]
        thresholds = list(itertools.product(pvalues, foldchanges))
        thresholds.pop(0)
    else:
        pvalues = ["1"]
        foldchanges = ["1.5","2"]
        thresholds = list(itertools.product(pvalues, foldchanges))

    for i in reversed(range(number_of_comparisons)):
        
        getFoldChange(i,history,history_client,dataDirectory,galaxy_host)

        
    for index, values in reversed(list(enumerate(list(itertools.product(range(number_of_comparisons),thresholds))))):
        
        (comparison,(pvalue, foldchange)) = values
                              
        print(index)
        print(values)
        
        
        getStringNetworks(index,comparison,history,history_client,dataDirectory,galaxy_host)
        getBioGridNetworks(index,comparison,history,history_client,dataDirectory,galaxy_host)
    
        
        getPathways(index,comparison,pvalue, foldchange,history,history_client,dataDirectory,galaxy_host)
 
        getDrugEnrichment(index,comparison,pvalue, foldchange,history,history_client,dataDirectory,galaxy_host)
        getGOEnrichment(index,comparison,pvalue, foldchange,history,history_client,dataDirectory,wwwDirectoryPlots, galaxy_host)
        
        if species in ["Human","Mouse"]:
            getTFs(index,comparison,pvalue, foldchange,history,history_client,dataDirectory,galaxy_host)

        


with open('data/expTable.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        getGalaxyData(row['ID'],row['Type'],row["Species"],row["foldChangeOnly"])