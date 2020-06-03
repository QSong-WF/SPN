
•	“MainCode.R”
The file is used for searching networks initiated from each gene located on the human functional protein interaction (FI) database. 
Input: the gene expression data of breast cancer patients; survival time of breast cancer patients; protein-protein interaction database;



•	“SignificanceTest.R”
The document is used for selecting significant survival prognostic network markers (SPNs) through the results of three 100 trials.
Input: the results of three 100 trials; the subnetworks obtained from “MainCode.R” initiated from each gene; corresponding scores of subnetworks;
Output: significant survival prognostic network markers (SPNs).


•	“SurvialAnalysis.R”
The document is used for survival analysis of patients based on identified subnetwork markes.
Input: SPNs identifed after the three tests; the gene expression data of breast cancer patients; survival time of breast cancer patients;
Output: the results of survival analysis including the survival curves, the classification of patients, and the P-value of log-rank test.

•	“FisherTest.R”
The document is the program used for Fisher exact test in BP sets and KEGG pathway sets.
Input: significant survival prognostic network markers (SPNs); KEGG pathway sets; Biological process (BP) sets;
Output: enrichment of SPNs in KEGG pathway sets and BP sets.


•	“MyImageFunction.R”
The document is used for plotting the enrichment results in BP sets and KEGG pathway sets.
Input: the enrichment result of SPNs in BP sets and KEGG sets;
Output: the plots showing the enrichment result of SPNs in BP sets and KEGG sets.

Datasets:
•	Training set
The breast cancer patients dataset used as training set can be downloaded from the additional file 1 in paper PMID: 23618380.

•	Test set
The breast cancer patients datasets used as test set can be downloaded from the GEO database, referenced by accession number GSE25066. 

•	“FIdatabase.txt”
The document of human functional protein interaction (FI) database can be downloaded from the paper PMID: 20482850.

•	“ProlferationGene.txt”
The document of proliferation genes used for classifying patients into different proliferation tertile can be downloaded from the paper PMID: 23618380.

•	“BP.txt” and “KEGG.txt”
The document of biological process sets and the KEGG pathway sets downloaded from the MSigDB database of GSEA website.
