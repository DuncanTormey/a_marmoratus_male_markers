
# _Aspidoscelis marmoratus_ Male Specific Markers

## Imports and Constants


```python
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib.patches import Rectangle
import seaborn as sns
import os
from glob import glob
from wand.image import Image as WImage
import itertools as it
from Bio import SeqIO
from sklearn.preprocessing import scale
from scipy.spatial import distance
from scipy.cluster import hierarchy
import primer3
import sys
import datetime

```


```python
sys.path = ['../../../python_tools/fasta_class/'] + sys.path
```


```python
%matplotlib inline
sns.set_style("whitegrid")
```


```python
# plotting constants
minor_f_size = 14
major_f_size = 16
title_f_size = 18
notebook_fig_size = (10,8)
sns.set_context(rc={"font.size":major_f_size,"axes.titlesize":major_f_size,"axes.labelsize":minor_f_size})   


# other constants
marm_genome_size = 1639530780
read_length = 250
```

## Gathering Data

This family of lizards was originally sequenced in order to generate a meotic map for the _A. marmoratus_ genome assembly. Since we already have the DNAseq data for the animals, we decided to try and leverage those data in order to identify a male specific dna sequence that can be used for genotyping. In this section I will download/link to the dnaseq alignments that Aaron Odell generated. However, first I downloaded some meta data as well as the data about the original sequencing run

### Sequencing Data

This cohort of animals was sequenced across three illumina HiSeq 2500 2 lane flowcells in rapid run mode.


```python
sampleReportPaths = glob('/n/analysis/Baumann/*/MOLNG-1807/*/Sample_Report.csv')
sampleReportPaths
```




    ['/n/analysis/Baumann/rmh/MOLNG-1807/H2L7VBCXYa/Sample_Report.csv',
     '/n/analysis/Baumann/rmh/MOLNG-1807/H7CN3BCXYa/Sample_Report.csv',
     '/n/analysis/Baumann/rmh/MOLNG-1807/H2WCCBCXYa/Sample_Report.csv']




```python
sampleReportDfs = [pd.read_csv(path) for path in sampleReportPaths]
sampleReportDf = pd.concat(sampleReportDfs)
sampleReportDf.columns = [c.replace(' ', '_') for c in sampleReportDf.columns] # I dont like to deal with spaces in column names
sampleReportDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>output</th>
      <th>order</th>
      <th>order_type</th>
      <th>lane</th>
      <th>sample_name</th>
      <th>library_id</th>
      <th>illumina_index</th>
      <th>custom_barcode</th>
      <th>read</th>
      <th>reference</th>
      <th>lab</th>
      <th>total_reads</th>
      <th>pass_filter_reads</th>
      <th>pass_filter_percent</th>
      <th>align_percent</th>
      <th>type</th>
      <th>read_length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>s_1_1_CGATGT.fastq.gz</td>
      <td>MOLNG-1807</td>
      <td>DNA-Seq (1ug)</td>
      <td>1</td>
      <td>A.tig_17115</td>
      <td>L24712</td>
      <td>CGATGT</td>
      <td>NaN</td>
      <td>1</td>
      <td>None</td>
      <td>Baumann Lab</td>
      <td>18234321</td>
      <td>18234321</td>
      <td>100.0</td>
      <td>NaN</td>
      <td>Paired Reads</td>
      <td>251</td>
    </tr>
    <tr>
      <th>1</th>
      <td>s_1_2_CGATGT.fastq.gz</td>
      <td>MOLNG-1807</td>
      <td>DNA-Seq (1ug)</td>
      <td>1</td>
      <td>A.tig_17115</td>
      <td>L24712</td>
      <td>CGATGT</td>
      <td>NaN</td>
      <td>2</td>
      <td>None</td>
      <td>Baumann Lab</td>
      <td>18234321</td>
      <td>18234321</td>
      <td>100.0</td>
      <td>NaN</td>
      <td>Paired Reads</td>
      <td>251</td>
    </tr>
    <tr>
      <th>2</th>
      <td>s_1_1_TGACCA.fastq.gz</td>
      <td>MOLNG-1807</td>
      <td>DNA-Seq (1ug)</td>
      <td>1</td>
      <td>A.tig_17118</td>
      <td>L24713</td>
      <td>TGACCA</td>
      <td>NaN</td>
      <td>1</td>
      <td>None</td>
      <td>Baumann Lab</td>
      <td>14017603</td>
      <td>14017603</td>
      <td>100.0</td>
      <td>NaN</td>
      <td>Paired Reads</td>
      <td>251</td>
    </tr>
    <tr>
      <th>3</th>
      <td>s_1_2_TGACCA.fastq.gz</td>
      <td>MOLNG-1807</td>
      <td>DNA-Seq (1ug)</td>
      <td>1</td>
      <td>A.tig_17118</td>
      <td>L24713</td>
      <td>TGACCA</td>
      <td>NaN</td>
      <td>2</td>
      <td>None</td>
      <td>Baumann Lab</td>
      <td>14017603</td>
      <td>14017603</td>
      <td>100.0</td>
      <td>NaN</td>
      <td>Paired Reads</td>
      <td>251</td>
    </tr>
    <tr>
      <th>4</th>
      <td>s_1_1_ACAGTG.fastq.gz</td>
      <td>MOLNG-1807</td>
      <td>DNA-Seq (1ug)</td>
      <td>1</td>
      <td>A.tig_21321</td>
      <td>L24714</td>
      <td>ACAGTG</td>
      <td>NaN</td>
      <td>1</td>
      <td>None</td>
      <td>Baumann Lab</td>
      <td>9493500</td>
      <td>9493500</td>
      <td>100.0</td>
      <td>NaN</td>
      <td>Paired Reads</td>
      <td>251</td>
    </tr>
  </tbody>
</table>
</div>




```python
indexToSampleDict = dict(zip(sampleReportDf.illumina_index, sampleReportDf.sample_name))
indexToSampleDict
```




    {'ACAGTG': 'A.tig_21321',
     'ACTTGA': 'A.tig_22151',
     'ATCACG': 'A.tig_21545',
     'CAGATC': 'A.tig_21323',
     'CGATGT': 'A.tig_17115',
     'CTTGTA': 'A.tig_21544',
     'GCCAAT': 'A.tig_21322',
     'TGACCA': 'A.tig_17118',
     'TTAGGC': 'A.tig_21546'}




```python
readCounts = sampleReportDf.groupby('sample_name').sum().reset_index()[['sample_name', 'total_reads']]
readCounts['total_read_pairs'] = readCounts['total_reads'] / 2
readCounts['estimated_coverage'] = (readCounts['total_reads'] * 250)/ marm_genome_size
readCounts
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_name</th>
      <th>total_reads</th>
      <th>total_read_pairs</th>
      <th>estimated_coverage</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>A.tig_17115</td>
      <td>233057978</td>
      <td>116528989.0</td>
      <td>35.537298</td>
    </tr>
    <tr>
      <th>1</th>
      <td>A.tig_17118</td>
      <td>179242302</td>
      <td>89621151.0</td>
      <td>27.331341</td>
    </tr>
    <tr>
      <th>2</th>
      <td>A.tig_21321</td>
      <td>124487662</td>
      <td>62243831.0</td>
      <td>18.982209</td>
    </tr>
    <tr>
      <th>3</th>
      <td>A.tig_21322</td>
      <td>152940748</td>
      <td>76470374.0</td>
      <td>23.320811</td>
    </tr>
    <tr>
      <th>4</th>
      <td>A.tig_21323</td>
      <td>193542766</td>
      <td>96771383.0</td>
      <td>29.511914</td>
    </tr>
    <tr>
      <th>5</th>
      <td>A.tig_21544</td>
      <td>210063366</td>
      <td>105031683.0</td>
      <td>32.031019</td>
    </tr>
    <tr>
      <th>6</th>
      <td>A.tig_21545</td>
      <td>186505874</td>
      <td>93252937.0</td>
      <td>28.438910</td>
    </tr>
    <tr>
      <th>7</th>
      <td>A.tig_21546</td>
      <td>287261598</td>
      <td>143630799.0</td>
      <td>43.802410</td>
    </tr>
    <tr>
      <th>8</th>
      <td>A.tig_22151</td>
      <td>241724858</td>
      <td>120862429.0</td>
      <td>36.858847</td>
    </tr>
  </tbody>
</table>
</div>




```python
ax = readCounts.plot('sample_name', 
                     'total_read_pairs', 
                     kind='bar', 
                     fontsize=minor_f_size, 
                     legend=False, 
                     figsize=notebook_fig_size)

ax.yaxis.offsetText.set_fontsize(minor_f_size)
ax.set_yticklabels(['{:,}'.format(int(y)) for y in ax.get_yticks().tolist()]);
ax.set_ylabel('Number of Read Pairs', fontsize=major_f_size)
ax.set_xlabel('')
ax.set_title('Read Pairs Sequenced For Each Animal', fontsize=title_f_size)
```




    <matplotlib.text.Text at 0x7f30cbf34470>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_14_1.png)



```python
ax = readCounts.plot('sample_name', 
                     'estimated_coverage', 
                     kind='bar', 
                     color='darkred', 
                     fontsize=minor_f_size, 
                     legend=False, 
                     figsize=notebook_fig_size)

ax.yaxis.offsetText.set_fontsize(minor_f_size)
ax.set_yticklabels(['{:,}'.format(int(y)) for y in ax.get_yticks().tolist()]);
ax.set_ylabel('Estimated Coverage', fontsize=major_f_size)
ax.set_xlabel('')
ax.set_title('Coverage', fontsize=title_f_size)
```




    <matplotlib.text.Text at 0x7f30cbe8aa20>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_15_1.png)


### Meta Data


```python
metaData = pd.read_excel('../data/external_data/mm_animal_meta_data_1.xlsx')
metaData.columns = [c.replace(' ', '_').lower() for c in metaData.columns]
metaData
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>animal_id</th>
      <th>sample_name</th>
      <th>generation</th>
      <th>letter</th>
      <th>species</th>
      <th>laid_by</th>
      <th>date_hatched</th>
      <th>hatchling_notes</th>
      <th>gender</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>21321</td>
      <td>A.tig_21321</td>
      <td>F1</td>
      <td>Z259</td>
      <td>tig</td>
      <td>17118.0</td>
      <td>2016-03-29</td>
      <td>NaN</td>
      <td>F</td>
    </tr>
    <tr>
      <th>1</th>
      <td>21322</td>
      <td>A.tig_21322</td>
      <td>F1</td>
      <td>Z259</td>
      <td>tig</td>
      <td>17118.0</td>
      <td>2016-03-29</td>
      <td>NaN</td>
      <td>M</td>
    </tr>
    <tr>
      <th>2</th>
      <td>21323</td>
      <td>A.tig_21323</td>
      <td>F1</td>
      <td>Z259</td>
      <td>tig</td>
      <td>17118.0</td>
      <td>2016-03-29</td>
      <td>NaN</td>
      <td>M</td>
    </tr>
    <tr>
      <th>3</th>
      <td>21544</td>
      <td>A.tig_21544</td>
      <td>F1</td>
      <td>C263</td>
      <td>tig</td>
      <td>17118.0</td>
      <td>2016-04-24</td>
      <td>NaN</td>
      <td>F</td>
    </tr>
    <tr>
      <th>4</th>
      <td>21545</td>
      <td>A.tig_21545</td>
      <td>F1</td>
      <td>C263</td>
      <td>tig</td>
      <td>17118.0</td>
      <td>2016-04-24</td>
      <td>NaN</td>
      <td>M</td>
    </tr>
    <tr>
      <th>5</th>
      <td>21546</td>
      <td>A.tig_21546</td>
      <td>F1</td>
      <td>C263</td>
      <td>tig</td>
      <td>17118.0</td>
      <td>2016-04-29</td>
      <td>NaN</td>
      <td>F</td>
    </tr>
    <tr>
      <th>6</th>
      <td>22151</td>
      <td>A.tig_22151</td>
      <td>F1</td>
      <td>I272</td>
      <td>tig</td>
      <td>17118.0</td>
      <td>2016-06-14</td>
      <td>hatched with opening in chest and kinked tail</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>7</th>
      <td>17118</td>
      <td>A.tig_17118</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>tig</td>
      <td>NaN</td>
      <td>NaT</td>
      <td>NaN</td>
      <td>F</td>
    </tr>
    <tr>
      <th>8</th>
      <td>17115</td>
      <td>A.tig_17115</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>tig</td>
      <td>NaN</td>
      <td>NaT</td>
      <td>NaN</td>
      <td>M</td>
    </tr>
  </tbody>
</table>
</div>




```python
sampleToSexDict = dict(zip(metaData.sample_name, metaData.gender))
```

### Alignment Files


```python
alignmentPaths = glob('/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/*_merged_sorted_markDup_indelRealign.bam')
alignmentPaths
```




    ['/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ATCACG_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACAGTG_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/GCCAAT_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CTTGTA_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CGATGT_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACTTGA_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TTAGGC_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TGACCA_merged_sorted_markDup_indelRealign.bam',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CAGATC_merged_sorted_markDup_indelRealign.bam']




```python
alignmentIndexPaths = glob('/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/*_merged_sorted_markDup_indelRealign.bai')
alignmentIndexPaths
```




    ['/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CGATGT_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TGACCA_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CAGATC_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CTTGTA_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACAGTG_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TTAGGC_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACTTGA_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ATCACG_merged_sorted_markDup_indelRealign.bai',
     '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/GCCAAT_merged_sorted_markDup_indelRealign.bai']



### Genome Assembly Reference Used


```python
genomeAssemblyPath = "/n/projects/aodell/GENOMES/a_marmorata_unmasked.fasta"
genomePicardDict = "/n/projects/aodell/GENOMES/a_marmorata_unmasked.fasta.fai"
```

## Alignment Mapping Stats

First I am going use the local install of samtools. Here is the info about the version I used.


```bash
%%bash
samtools help | head
```

    
    Program: samtools (Tools for alignments in the SAM format)
    Version: 1.3.1 (using htslib 1.3.1)
    
    Usage:   samtools <command> [options]
    
    Commands:
      -- Indexing
         dict           create a sequence dictionary file
         faidx          index/extract FASTA


 Write the commands to a text file so I can parallelize the execution.


```python
flagstatPathsDict = {}
with open('../bin/samtools_flagstat.txt','w') as fw:
    for bam in alignmentPaths:
        filename = os.path.basename(bam)
        flagstat = '../data/' + filename.replace('.bam', '.flagstat')
        index = filename.split("_")[0]
        sample_name = indexToSampleDict[index]
        command = 'samtools flagstat {} > {}'.format(bam, flagstat)
        flagstatPathsDict[sample_name] = flagstat
        fw.write(command + '\n')

print(len(flagstatPathsDict.keys()))
```

    9


Run the commands in parallel. Took about 5 min. (I actually ran the command in terminal).


```bash
%%bash
cat ../bin/samtools_flagstat.txt
#nohup parallel --no-notice -j 9 :::: ../bin/samtools_flagstat.txt &> ../bin/samtools_flagstat.out &
```

    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ATCACG_merged_sorted_markDup_indelRealign.bam > ../data/ATCACG_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACAGTG_merged_sorted_markDup_indelRealign.bam > ../data/ACAGTG_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/GCCAAT_merged_sorted_markDup_indelRealign.bam > ../data/GCCAAT_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CTTGTA_merged_sorted_markDup_indelRealign.bam > ../data/CTTGTA_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CGATGT_merged_sorted_markDup_indelRealign.bam > ../data/CGATGT_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACTTGA_merged_sorted_markDup_indelRealign.bam > ../data/ACTTGA_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TTAGGC_merged_sorted_markDup_indelRealign.bam > ../data/TTAGGC_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TGACCA_merged_sorted_markDup_indelRealign.bam > ../data/TGACCA_merged_sorted_markDup_indelRealign.flagstat
    samtools flagstat /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CAGATC_merged_sorted_markDup_indelRealign.bam > ../data/CAGATC_merged_sorted_markDup_indelRealign.flagstat


Clean up the results file so that I can read them in as a csv


```bash
%%bash 
rm ../data/flagstat_mapped_data.csv
touch ../data/flagstat_mapped_data.csv
for file in ../data/*Realign.flagstat
 do
  echo $file
  index=$(basename $file | cut -f 1 -d '_')
  cat $file | grep "mapped (" | sed "s/ + 0 mapped (/,/" | sed "s/%.*$//" | sed "s/^/$index,/" >> ../data/flagstat_mapped_data.csv
 done

head -n 5 ../data/flagstat_mapped_data.csv
```

    ../data/ACAGTG_merged_sorted_markDup_indelRealign.flagstat
    ../data/ACTTGA_merged_sorted_markDup_indelRealign.flagstat
    ../data/ATCACG_merged_sorted_markDup_indelRealign.flagstat
    ../data/CAGATC_merged_sorted_markDup_indelRealign.flagstat
    ../data/CGATGT_merged_sorted_markDup_indelRealign.flagstat
    ../data/CTTGTA_merged_sorted_markDup_indelRealign.flagstat
    ../data/GCCAAT_merged_sorted_markDup_indelRealign.flagstat
    ../data/TGACCA_merged_sorted_markDup_indelRealign.flagstat
    ../data/TTAGGC_merged_sorted_markDup_indelRealign.flagstat
    ACAGTG,121793509,94.99
    ACTTGA,238250452,97.18
    ATCACG,185085996,96.62
    CAGATC,193085899,97.05
    CGATGT,228954817,96.85



```python
flagstatDF = pd.read_csv('../data/flagstat_mapped_data.csv',names=['index', 'reads_mapped', 'percent_mapped'])
flagstatDF['sample_name'] = flagstatDF['index'].apply(lambda x: indexToSampleDict[x])
flagstatDF['sex'] = flagstatDF['sample_name'].apply(lambda x: sampleToSexDict[x])
flagstatDF.sort_values('sex', inplace=True)
flagstatDF
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>index</th>
      <th>reads_mapped</th>
      <th>percent_mapped</th>
      <th>sample_name</th>
      <th>sex</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ACAGTG</td>
      <td>121793509</td>
      <td>94.99</td>
      <td>A.tig_21321</td>
      <td>F</td>
    </tr>
    <tr>
      <th>5</th>
      <td>CTTGTA</td>
      <td>205882982</td>
      <td>96.72</td>
      <td>A.tig_21544</td>
      <td>F</td>
    </tr>
    <tr>
      <th>7</th>
      <td>TGACCA</td>
      <td>179237242</td>
      <td>97.17</td>
      <td>A.tig_17118</td>
      <td>F</td>
    </tr>
    <tr>
      <th>8</th>
      <td>TTAGGC</td>
      <td>282106059</td>
      <td>97.50</td>
      <td>A.tig_21546</td>
      <td>F</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ATCACG</td>
      <td>185085996</td>
      <td>96.62</td>
      <td>A.tig_21545</td>
      <td>M</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CAGATC</td>
      <td>193085899</td>
      <td>97.05</td>
      <td>A.tig_21323</td>
      <td>M</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CGATGT</td>
      <td>228954817</td>
      <td>96.85</td>
      <td>A.tig_17115</td>
      <td>M</td>
    </tr>
    <tr>
      <th>6</th>
      <td>GCCAAT</td>
      <td>150815423</td>
      <td>96.03</td>
      <td>A.tig_21322</td>
      <td>M</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ACTTGA</td>
      <td>238250452</td>
      <td>97.18</td>
      <td>A.tig_22151</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
df = flagstatDF.dropna().copy()
df['color'] = df['sex'].apply(lambda x: 'purple' if x=='M' else 'green')
ax = df.plot('sample_name', 
             'percent_mapped', 
             kind='bar', 
             color=df['color'], 
             fontsize=minor_f_size,
             legend=False, 
             figsize=notebook_fig_size)

ax.set_ylim(94,98)
ax.yaxis.offsetText.set_fontsize(minor_f_size)
ax.set_yticklabels(['{:,}'.format(int(y)) for y in ax.get_yticks().tolist()]);
ax.set_ylabel('Percent Mapped', fontsize=major_f_size)
ax.set_xlabel('')


#legend
female = plt.Rectangle((0,0),1,1,fc="purple", edgecolor = 'none')
male = plt.Rectangle((0,0),1,1,fc='green',  edgecolor = 'none')

l = ax.legend([female, male], 
              ['F', 'M'], 
              loc=1, 
              ncol = 1, 
              fontsize=minor_f_size
              )

l.draw_frame(False)
```


![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_33_0.png)


## Extracting Unmapped Reads


```python
unmappedPathsDict = {}
with open('../bin/extract_unmapped.txt','w') as fw:
    for bam in alignmentPaths:
        filename = os.path.basename(bam)
        unmapped = '../data/' + filename.replace('.bam', '.unmapped.bam')
        index = filename.split("_")[0]
        sample_name = indexToSampleDict[index]
        command = 'samtools view -@ 3 -f 12 -b {} > {}'.format(bam, unmapped)
        unmappedPathsDict[sample_name] = unmapped
        fw.write(command + '\n')

print(len(unmappedPathsDict.keys()))
```

    9



```bash
%%bash
cat ../bin/extract_unmapped.txt
#nohup parallel --no-notice -j 27 :::: ../bin/extract_unmapped.txt &> ../bin/extract_unmapped.out &
```

    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ATCACG_merged_sorted_markDup_indelRealign.bam > ../data/ATCACG_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACAGTG_merged_sorted_markDup_indelRealign.bam > ../data/ACAGTG_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/GCCAAT_merged_sorted_markDup_indelRealign.bam > ../data/GCCAAT_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CTTGTA_merged_sorted_markDup_indelRealign.bam > ../data/CTTGTA_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CGATGT_merged_sorted_markDup_indelRealign.bam > ../data/CGATGT_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACTTGA_merged_sorted_markDup_indelRealign.bam > ../data/ACTTGA_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TTAGGC_merged_sorted_markDup_indelRealign.bam > ../data/TTAGGC_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TGACCA_merged_sorted_markDup_indelRealign.bam > ../data/TGACCA_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools view -@ 3 -f 12 -b /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CAGATC_merged_sorted_markDup_indelRealign.bam > ../data/CAGATC_merged_sorted_markDup_indelRealign.unmapped.bam



```python
sexPathsDict = {}
with open('../bin/merge_sexes.txt','w') as fw:
    male = []
    female = []
    for sample_name in unmappedPathsDict:
        if sampleToSexDict[sample_name]=='M':
            male.append(unmappedPathsDict[sample_name])
        elif sampleToSexDict[sample_name]=='F':
            female.append(unmappedPathsDict[sample_name])
        else:
            pass

            
    male_bam = '../data/male_unmapped.bam'
    female_bam = '../data/female_unmapped.bam'
    command = 'samtools merge -f {} {}'.format(male_bam, ' '.join(male))
    command2 = 'samtools merge -f {} {}'.format(female_bam, ' '.join(female))
    sexPathsDict['M'] = female_bam
    sexPathsDict['F'] = male_bam
    
    fw.write(command + '\n' + command2 + '\n')


```


```bash
%%bash
cat ../bin/merge_sexes.txt
#nohup parallel --no-notice -j 2 :::: ../bin/merge_sexes.txt &> ../bin/merge_sexes.out &
```

    samtools merge -f ../data/male_unmapped.bam ../data/ATCACG_merged_sorted_markDup_indelRealign.unmapped.bam ../data/GCCAAT_merged_sorted_markDup_indelRealign.unmapped.bam ../data/CGATGT_merged_sorted_markDup_indelRealign.unmapped.bam ../data/CAGATC_merged_sorted_markDup_indelRealign.unmapped.bam
    samtools merge -f ../data/female_unmapped.bam ../data/ACAGTG_merged_sorted_markDup_indelRealign.unmapped.bam ../data/CTTGTA_merged_sorted_markDup_indelRealign.unmapped.bam ../data/TTAGGC_merged_sorted_markDup_indelRealign.unmapped.bam ../data/TGACCA_merged_sorted_markDup_indelRealign.unmapped.bam


## Alignments Stats for Assemblies

Download Picard Tools for insert statistics.


```bash
%%bash
cd ../bin
#git clone https://github.com/broadinstitute/picard.git
#cd picard/
#./gradlew shadowJar
#java -jar build/libs/picard.jar
```

Sub sample alignments to gather insert size metrics


```python
alignmentPathsDict = {indexToSampleDict[os.path.basename(path).split('_')[0]]:path for path in alignmentPaths}
alignmentPathsDict
```




    {'A.tig_17115': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CGATGT_merged_sorted_markDup_indelRealign.bam',
     'A.tig_17118': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TGACCA_merged_sorted_markDup_indelRealign.bam',
     'A.tig_21321': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACAGTG_merged_sorted_markDup_indelRealign.bam',
     'A.tig_21322': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/GCCAAT_merged_sorted_markDup_indelRealign.bam',
     'A.tig_21323': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CAGATC_merged_sorted_markDup_indelRealign.bam',
     'A.tig_21544': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CTTGTA_merged_sorted_markDup_indelRealign.bam',
     'A.tig_21545': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ATCACG_merged_sorted_markDup_indelRealign.bam',
     'A.tig_21546': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TTAGGC_merged_sorted_markDup_indelRealign.bam',
     'A.tig_22151': '/n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACTTGA_merged_sorted_markDup_indelRealign.bam'}




```python
sexSubPathsDict = {'F':[], 'M':[]}

with open('../bin/sub_all_sexes.txt','w') as fw:
    for sample_name in alignmentPathsDict:
        if sampleToSexDict[sample_name]=='M' or sampleToSexDict[sample_name]=='F':
            mapped_reads = flagstatDF[flagstatDF['sample_name'] == sample_name]['reads_mapped'].iloc[0]
            perc = round(250e3/mapped_reads,6)
            out_bam = '../data/' + os.path.basename(alignmentPathsDict[sample_name]).replace('.bam','.250kSub.bam')
            command = 'samtools view -@ 2 -b -F 4 -s {} {} > {}'.format(perc, alignmentPathsDict[sample_name], out_bam)
            print(command)
            
            fw.write(command + '\n')
            sexSubPathsDict[sampleToSexDict[sample_name]].append(out_bam)
            
sexSubPathsDict   
```

    samtools view -@ 2 -b -F 4 -s 0.001351 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ATCACG_merged_sorted_markDup_indelRealign.bam > ../data/ATCACG_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.002053 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACAGTG_merged_sorted_markDup_indelRealign.bam > ../data/ACAGTG_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001658 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/GCCAAT_merged_sorted_markDup_indelRealign.bam > ../data/GCCAAT_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001214 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CTTGTA_merged_sorted_markDup_indelRealign.bam > ../data/CTTGTA_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001092 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CGATGT_merged_sorted_markDup_indelRealign.bam > ../data/CGATGT_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.000886 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TTAGGC_merged_sorted_markDup_indelRealign.bam > ../data/TTAGGC_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001395 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TGACCA_merged_sorted_markDup_indelRealign.bam > ../data/TGACCA_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001295 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CAGATC_merged_sorted_markDup_indelRealign.bam > ../data/CAGATC_merged_sorted_markDup_indelRealign.250kSub.bam





    {'F': ['../data/ACAGTG_merged_sorted_markDup_indelRealign.250kSub.bam',
      '../data/CTTGTA_merged_sorted_markDup_indelRealign.250kSub.bam',
      '../data/TTAGGC_merged_sorted_markDup_indelRealign.250kSub.bam',
      '../data/TGACCA_merged_sorted_markDup_indelRealign.250kSub.bam'],
     'M': ['../data/ATCACG_merged_sorted_markDup_indelRealign.250kSub.bam',
      '../data/GCCAAT_merged_sorted_markDup_indelRealign.250kSub.bam',
      '../data/CGATGT_merged_sorted_markDup_indelRealign.250kSub.bam',
      '../data/CAGATC_merged_sorted_markDup_indelRealign.250kSub.bam']}




```bash
%%bash
cat ../bin/sub_all_sexes.txt
#nohup parallel --no-notice -j 16 :::: ../bin/sub_all_sexes.txt &> ../bin/sub_all_sexes.out &
```

    samtools view -@ 2 -b -F 4 -s 0.001351 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ATCACG_merged_sorted_markDup_indelRealign.bam > ../data/ATCACG_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.002053 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/ACAGTG_merged_sorted_markDup_indelRealign.bam > ../data/ACAGTG_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001658 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/GCCAAT_merged_sorted_markDup_indelRealign.bam > ../data/GCCAAT_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001214 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CTTGTA_merged_sorted_markDup_indelRealign.bam > ../data/CTTGTA_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001092 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CGATGT_merged_sorted_markDup_indelRealign.bam > ../data/CGATGT_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.000886 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TTAGGC_merged_sorted_markDup_indelRealign.bam > ../data/TTAGGC_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001395 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/TGACCA_merged_sorted_markDup_indelRealign.bam > ../data/TGACCA_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools view -@ 2 -b -F 4 -s 0.001295 /n/projects/aodell/A_MARMORATA_MEIOTIC_MAP/BWA/UNIFIED_GENOTYPER_BAMS/CAGATC_merged_sorted_markDup_indelRealign.bam > ../data/CAGATC_merged_sorted_markDup_indelRealign.250kSub.bam


Merge the subsampled alignments.


```python
subSexMergeCommand1 = 'samtools merge -@ 8 -f ../data/subSexMale.bam {}'.format(' '.join(sexSubPathsDict['M']))
subSexMergeCommand2 = 'samtools merge -@ 8 -f ../data/subSexFemale.bam {}'.format(' '.join(sexSubPathsDict['F']))
```


```bash
%%bash -s "$subSexMergeCommand1" "$subSexMergeCommand2 "
echo $1
echo $2
#eval $1
#eval $2
```

    samtools merge -@ 8 -f ../data/subSexMale.bam ../data/ATCACG_merged_sorted_markDup_indelRealign.250kSub.bam ../data/GCCAAT_merged_sorted_markDup_indelRealign.250kSub.bam ../data/CGATGT_merged_sorted_markDup_indelRealign.250kSub.bam ../data/CAGATC_merged_sorted_markDup_indelRealign.250kSub.bam
    samtools merge -@ 8 -f ../data/subSexFemale.bam ../data/ACAGTG_merged_sorted_markDup_indelRealign.250kSub.bam ../data/CTTGTA_merged_sorted_markDup_indelRealign.250kSub.bam ../data/TTAGGC_merged_sorted_markDup_indelRealign.250kSub.bam ../data/TGACCA_merged_sorted_markDup_indelRealign.250kSub.bam


Run Picard Tools on the merged subsampled alignments.


```bash
%%bash
#java -jar ../bin/picard/build/libs/picard.jar CollectInsertSizeMetrics I=../data/subSexMale.bam O=../data/subSexMale.insert_metrics.txt H=../data/subSexMale.insert_metrics.hist.pdf M=0.05
```


```bash
%%bash
#java -jar ../bin/picard/build/libs/picard.jar CollectInsertSizeMetrics I=../data/../data/subSexFemale.bam O=../data/subSexFemale.insert_metrics.txt H=../data/subSexFemale.insert_metrics.hist.pdf M=0.05

```

Plots of insert size and metrics output by picard tools.


```python
WImage(filename='../data/../data/subSexMale.insert_metrics.hist.pdf')
```




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_53_0.png)




```python
WImage(filename='../data/subSexFemale.insert_metrics.hist.pdf')
```




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_54_0.png)




```bash
%%bash
grep -A 2 "## METRICS CLASS" ../data/subSexFemale.insert_metrics.txt | sed 1d | cut -f 1-6| column -t
```

    MEDIAN_INSERT_SIZE  MEDIAN_ABSOLUTE_DEVIATION  MIN_INSERT_SIZE  MAX_INSERT_SIZE  MEAN_INSERT_SIZE  STANDARD_DEVIATION
    535                 31                         2                72207764         523.151352        86.359569



```bash
%%bash
grep -A 2 "## METRICS CLASS" ../data/subSexMale.insert_metrics.txt | sed 1d | cut -f 1-6| column -t
```

    MEDIAN_INSERT_SIZE  MEDIAN_ABSOLUTE_DEVIATION  MIN_INSERT_SIZE  MAX_INSERT_SIZE  MEAN_INSERT_SIZE  STANDARD_DEVIATION
    538                 30                         2                72335463         531.725321        72.614146



```python
femaleInsertSize = 535
maleInsertSize = 538
```

## Assemblying Unmapped Reads 127-mer


```bash
%%bash
cd ../bin
wget https://sourceforge.net/projects/soapdenovo2/files/SOAPdenovo2/bin/r240/SOAPdenovo2-bin-LINUX-generic-r240.tgz
```


```bash
%%bash
cd ../bin
tar -zxf SOAPdenovo2-bin-LINUX-generic-r240.tgz
```

#### All Male Unmapped Reads Assembly 127-mer

First I need to make the config file for the SOAP run. I copied the text from the manual and modified it in the next cell.


```bash
%%bash
cd ../data
ls -1 $PWD/*.bam
```

    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ACAGTG_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ACAGTG_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ACTTGA_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ATCACG_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ATCACG_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_17115_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_17118_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21321_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21322_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21323_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21544_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21545_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21546_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CAGATC_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CAGATC_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CGATGT_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CGATGT_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CTTGTA_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CTTGTA_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/female_unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/F_rF.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/GCCAAT_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/GCCAAT_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/male_unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/M_rF.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/subSexFemale.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/subSexMale.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TGACCA_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TGACCA_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TTAGGC_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TTAGGC_merged_sorted_markDup_indelRealign.unmapped.bam



```bash
%%bash
#mkdir ../data/malelib
cd ../data/malelib
#ln -s ../male_unmapped.bam
ls $PWD/*
```

    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/malelib/male_unmapped.bam



```python
maleConfig ="""
maximal read length
max_rd_len=250
[LIB]
#average insert size
avg_ins=538
#if sequence needs to be reversed, 1 for forward reverse
reverse_seq=1
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=250
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#bam file for single or paired reads, reads 1 in paired reads file should always be followed by reads 2
#NOTE: If a read in bam file fails platform/vendor quality checks(the flag field 0x0200 is set), itself and it's paired read would be ignored.
b=/n/projects/dut/a_marmorata/dnaseq_sex_determ/data/malelib/male_unmapped.bam
"""
with open('../bin/male_Assembly.config','w') as fw:
    fw.write(maleConfig)
print(maleConfig)
```

    
    maximal read length
    max_rd_len=250
    [LIB]
    #average insert size
    avg_ins=538
    #if sequence needs to be reversed, 1 for forward reverse
    reverse_seq=1
    #in which part(s) the reads are used
    asm_flags=3
    #use only first 100 bps of each read
    rd_len_cutoff=250
    #in which order the reads are used while scaffolding
    rank=1
    # cutoff of pair number for a reliable connection (at least 3 for short insert size)
    pair_num_cutoff=3
    #minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
    map_len=32
    #bam file for single or paired reads, reads 1 in paired reads file should always be followed by reads 2
    #NOTE: If a read in bam file fails platform/vendor quality checks(the flag field 0x0200 is set), itself and it's paired read would be ignored.
    b=/n/projects/dut/a_marmorata/dnaseq_sex_determ/data/malelib/male_unmapped.bam
    



```bash
%%bash
cd ../bin/SOAPdenovo2-bin-LINUX-generic-r240
#nohup ./SOAPdenovo-127mer all -s ../male_Assembly.config -K 127 -R -p 20 -o male_assembly 1>ass.log 2>ass.err &
## after the assembly finished...
#mkdir ../../data/male_assembly
#mv male_assembly* ../../data/male_assembly/
#mv ass* ../../data/male_assembly/
```

##### All Male Kmer Freq


```python
maleKmerFreqDF= pd.read_csv('../data/male_assembly/male_assembly.kmerFreq', sep ='\t', names=['Number of Distinct Kmers'])
maleKmerFreqDF = maleKmerFreqDF.reset_index()
maleKmerFreqDF.columns = ['Frequency', 'Number of Distinct Kmers']
maleKmerFreqDF['Frequency'] = maleKmerFreqDF['Frequency'] + 1 
maleKmerFreqDF.head()
maleKmerFreqDF.tail()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Frequency</th>
      <th>Number of Distinct Kmers</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>250</th>
      <td>251</td>
      <td>1</td>
    </tr>
    <tr>
      <th>251</th>
      <td>252</td>
      <td>111</td>
    </tr>
    <tr>
      <th>252</th>
      <td>253</td>
      <td>0</td>
    </tr>
    <tr>
      <th>253</th>
      <td>254</td>
      <td>0</td>
    </tr>
    <tr>
      <th>254</th>
      <td>255</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
ax = maleKmerFreqDF[(maleKmerFreqDF['Frequency'] > 2) & (maleKmerFreqDF['Frequency'] < 100)].plot(x='Frequency', 
                                                                                                  y='Number of Distinct Kmers', 
                                                                                                  kind='bar',
                                                                                                  figsize=notebook_fig_size,
                                                                                                  fontsize=major_f_size,
                                                                                                  legend=False,
                                                                                                  title='Male 127-mer Assembly',
                                                                                                  color='green')
ax.set_xticks(np.arange(3,102,3));
ax.set_ylabel('Number of Distinct Kmers')
```




    <matplotlib.text.Text at 0x7f30cbcf1ef0>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_69_1.png)


##### Scaffold Lengths


```python
def ret_fasta_len(fasta_file):
    len_list = []
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        data = {'scaffold':rec.id,
               'length':len(rec.seq)}
        len_list.append(data)
    len_df = pd.DataFrame(len_list)
    len_df.sort_values('length')
    return len_df
```


```python
maleScaffoldLensDF = ret_fasta_len('../data/male_assembly/male_assembly.scafSeq')
ax = sns.distplot(maleScaffoldLensDF.length, kde=False)
ax.set_title('Male Assembly Scaffold Sizes')
```




    <matplotlib.text.Text at 0x7f30c73367b8>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_72_1.png)



```python
maleScaffoldLensDF.sort_values('length', ascending=False).head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>length</th>
      <th>scaffold</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>33</th>
      <td>4646</td>
      <td>scaffold34</td>
    </tr>
    <tr>
      <th>35</th>
      <td>3760</td>
      <td>scaffold36</td>
    </tr>
    <tr>
      <th>3445</th>
      <td>2726</td>
      <td>C7127</td>
    </tr>
    <tr>
      <th>3444</th>
      <td>2525</td>
      <td>C7125</td>
    </tr>
    <tr>
      <th>3443</th>
      <td>2105</td>
      <td>C7123</td>
    </tr>
  </tbody>
</table>
</div>



#### All Female Unmapped Reads Assembly

First I need to make the config file for the SOAP run. I copied the text from the manual and modified it in the next cell.


```bash
%%bash
cd ../data
ls -1 $PWD/*.bam
```

    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ACAGTG_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ACAGTG_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ACTTGA_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ATCACG_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/ATCACG_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_17115_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_17118_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21321_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21322_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21323_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21544_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21545_M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/A.tig_21546_F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CAGATC_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CAGATC_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CGATGT_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CGATGT_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CTTGTA_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/CTTGTA_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/female_unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/F_rF.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/F_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/GCCAAT_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/GCCAAT_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/male_unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/M_rF.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/M_rM.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/subSexFemale.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/subSexMale.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TGACCA_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TGACCA_merged_sorted_markDup_indelRealign.unmapped.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TTAGGC_merged_sorted_markDup_indelRealign.250kSub.bam
    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/TTAGGC_merged_sorted_markDup_indelRealign.unmapped.bam



```bash
%%bash
#mkdir ../data/femalelib
cd ../data/femalelib
#ln -s ../female_unmapped.bam
ls $PWD/*
```

    /n/projects/dut/a_marmorata/dnaseq_sex_determ/data/femalelib/female_unmapped.bam



```python
femaleConfig ="""
maximal read length
max_rd_len=250
[LIB]
#average insert size
avg_ins=538
#if sequence needs to be reversed, 1 for forward reverse
reverse_seq=1
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=250
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#bam file for single or paired reads, reads 1 in paired reads file should always be followed by reads 2
#NOTE: If a read in bam file fails platform/vendor quality checks(the flag field 0x0200 is set), itself and it's paired read would be ignored.
b=/n/projects/dut/a_marmorata/dnaseq_sex_determ/data/femalelib/female_unmapped.bam
"""
with open('../bin/female_Assembly.config','w') as fw:
    fw.write(femaleConfig)
print(femaleConfig)
```

    
    maximal read length
    max_rd_len=250
    [LIB]
    #average insert size
    avg_ins=538
    #if sequence needs to be reversed, 1 for forward reverse
    reverse_seq=1
    #in which part(s) the reads are used
    asm_flags=3
    #use only first 100 bps of each read
    rd_len_cutoff=250
    #in which order the reads are used while scaffolding
    rank=1
    # cutoff of pair number for a reliable connection (at least 3 for short insert size)
    pair_num_cutoff=3
    #minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
    map_len=32
    #bam file for single or paired reads, reads 1 in paired reads file should always be followed by reads 2
    #NOTE: If a read in bam file fails platform/vendor quality checks(the flag field 0x0200 is set), itself and it's paired read would be ignored.
    b=/n/projects/dut/a_marmorata/dnaseq_sex_determ/data/femalelib/female_unmapped.bam
    



```bash
%%bash
cd ../bin/SOAPdenovo2-bin-LINUX-generic-r240
#nohup ./SOAPdenovo-127mer all -s ../female_Assembly.config -K 127 -R -p 20 -o female_assembly 1>ass.log 2>ass.err &
## after the assembly finished...
#mkdir ../../data/female_assembly
#mv female_assembly* ../../data/female_assembly/
#mv ass* ../../data/female_assembly/
```

##### All Female Kmer Freq


```python
femaleKmerFreqDF= pd.read_csv('../data/female_assembly/female_assembly.kmerFreq', sep ='\t', names=['Number of Distinct Kmers'])
femaleKmerFreqDF = femaleKmerFreqDF.reset_index()
femaleKmerFreqDF.columns = ['Frequency', 'Number of Distinct Kmers']
femaleKmerFreqDF['Frequency'] = femaleKmerFreqDF['Frequency'] + 1 
femaleKmerFreqDF.head()
femaleKmerFreqDF.tail()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Frequency</th>
      <th>Number of Distinct Kmers</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>250</th>
      <td>251</td>
      <td>98</td>
    </tr>
    <tr>
      <th>251</th>
      <td>252</td>
      <td>7314</td>
    </tr>
    <tr>
      <th>252</th>
      <td>253</td>
      <td>0</td>
    </tr>
    <tr>
      <th>253</th>
      <td>254</td>
      <td>0</td>
    </tr>
    <tr>
      <th>254</th>
      <td>255</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
ax = femaleKmerFreqDF[(femaleKmerFreqDF['Frequency'] > 2) & (femaleKmerFreqDF['Frequency'] < 100)].plot(x='Frequency', 
                                                                                                        y='Number of Distinct Kmers', 
                                                                                                        kind='bar',
                                                                                                        figsize=notebook_fig_size,
                                                                                                        fontsize=major_f_size,
                                                                                                        legend=False,
                                                                                                        title='Female 127-mer Assembly',
                                                                                                        color='purple')
ax.set_xticks(np.arange(3,102,3));
ax.set_ylabel('Number of Distinct Kmers')
```




    <matplotlib.text.Text at 0x7f30c71f4d68>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_82_1.png)


##### Scaffold Sizes


```python
femaleScaffoldLensDF = ret_fasta_len('../data/female_assembly/female_assembly.scafSeq')
ax = sns.distplot(femaleScaffoldLensDF.length, kde=False)
ax.set_title('Female Assembly Scaffold Sizes')
ax.legend()
```

    /home/dut/anaconda3/lib/python3.6/site-packages/matplotlib/axes/_axes.py:545: UserWarning: No labelled objects found. Use label='...' kwarg on individual plots.
      warnings.warn("No labelled objects found. "



![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_84_1.png)



```python
femaleScaffoldLensDF.length.sum(), maleScaffoldLensDF.length.sum()
```




    (1901959, 1774579)




```python
ax = sns.distplot(femaleScaffoldLensDF.length, kde=False, label='Female')
ax = sns.distplot(maleScaffoldLensDF.length, kde=False, ax=ax, label='Male')
ax.set_title('Female and Male Assembly Scaffold Sizes')
ax.legend()
```




    <matplotlib.legend.Legend at 0x7f30c6ce7f60>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_86_1.png)


### Aligning All Male and All Female Reads to Assemblies

I will now align unmapped reads back to the assemblies. To do this I will use the system installation of bwa mem.


```bash
%%bash
bwa
```

    
    Program: bwa (alignment via Burrows-Wheeler transformation)
    Version: 0.7.15-r1140
    Contact: Heng Li <lh3@sanger.ac.uk>
    
    Usage:   bwa <command> [options]
    
    Command: index         index sequences in the FASTA format
             mem           BWA-MEM algorithm
             fastmap       identify super-maximal exact matches
             pemerge       merge overlapping paired ends (EXPERIMENTAL)
             aln           gapped/ungapped alignment
             samse         generate alignment (single ended)
             sampe         generate alignment (paired ended)
             bwasw         BWA-SW for long queries
    
             shm           manage indices in shared memory
             fa2pac        convert FASTA to PAC format
             pac2bwt       generate BWT from PAC
             pac2bwtgen    alternative algorithm for generating BWT
             bwtupdate     update .bwt to the new format
             bwt2sa        generate SA from BWT and Occ
    
    Note: To use BWA, you need to first index the genome with `bwa index'.
          There are three alignment algorithms in BWA: `mem', `bwasw', and
          `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
          first. Please `man ./bwa.1' for the manual.
    


Because bwa only takes fastq files as input, I will need to extract the reads from my unmapped bams. To do this I will use bedtools.


```bash
%%bash
bedtools bamtofastq
```

    
    Tool:    bedtools bamtofastq (aka bamToFastq)
    Version: v2.26.0
    Summary: Convert BAM alignments to FASTQ files. 
    
    Usage:   bamToFastq [OPTIONS] -i <BAM> -fq <FQ> 
    
    Options:
    	-fq2	FASTQ for second end.  Used if BAM contains paired-end data.
    		BAM should be sorted by query name is creating paired FASTQ.
    
    	-tags	Create FASTQ based on the mate info
    		in the BAM R2 and Q2 tags.
    
    Tips: 
    	If you want to create a single, interleaved FASTQ file 
    	for paired-end data, you can just write both to /dev/stdout:
    
    	bedtools bamtofastq -i x.bam -fq /dev/stdout -fq2 /dev/stdout > x.ilv.fq
    


####  Converting to FQ


```python
sexPathsDict
```




    {'F': '../data/male_unmapped.bam', 'M': '../data/female_unmapped.bam'}




```bash
%%bash
#bedtools bamtofastq -i ../data/male_unmapped.bam -fq ../data/male_unmapped_1.fq -fq2 ../data/male_unmapped_2.fq
```


```bash
%%bash
#bedtools bamtofastq -i ../data/female_unmapped.bam -fq ../data/female_unmapped_1.fq -fq2 ../data/female_unmapped_2.fq
```


```python
sexFqPaths = {
    'F':['../data/female_unmapped_1.fq', '../data/female_unmapped_2.fq'],
    'M':['../data/male_unmapped_1.fq', '../data/male_unmapped_2.fq']
             }
```

#### Index Reference


```python
assemblyPaths = {
    'rF':'../data/female_assembly/female_assembly.scafSeq',
    'rM':'../data/male_assembly/male_assembly.scafSeq'
}
```


```bash
%%bash
cd ../data/female_assembly/
#bwa index female_assembly.scafSeq
```


```bash
%%bash
cd ../data/male_assembly/
#bwa index male_assembly.scafSeq
```

#### Aligning 


```python
alignCombos = [(f, r) for f in sexFqPaths.keys() for r in assemblyPaths.keys()]
compBams = {}
with open('../bin/m_to_f_align.txt', 'w') as fw:
    for combo in alignCombos:
        out_bam = '../data/{}.bam'.format('_'.join(combo))
        temp = '../data/temp_{}'.format('_'.join(combo))
        bwa_command = 'bwa mem -t 8 {} {} {} | samtools sort -T {} -o {}'.format(
            assemblyPaths[combo[1]],
            sexFqPaths[combo[0]][0],
            sexFqPaths[combo[0]][1],
            temp,
            out_bam
        )
        compBams[combo] = out_bam
        fw.write(bwa_command + '\n')
```


```bash
%%bash
cat ../bin/m_to_f_align.txt
#parallel --no-notice -j 2 :::: ../bin/m_to_f_align.txt &> ../bin/m_to_f_align.out
```

    bwa mem -t 8 ../data/female_assembly/female_assembly.scafSeq ../data/female_unmapped_1.fq ../data/female_unmapped_2.fq | samtools sort -T ../data/temp_F_rF -o ../data/F_rF.bam
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/female_unmapped_1.fq ../data/female_unmapped_2.fq | samtools sort -T ../data/temp_F_rM -o ../data/F_rM.bam
    bwa mem -t 8 ../data/female_assembly/female_assembly.scafSeq ../data/male_unmapped_1.fq ../data/male_unmapped_2.fq | samtools sort -T ../data/temp_M_rF -o ../data/M_rF.bam
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/male_unmapped_1.fq ../data/male_unmapped_2.fq | samtools sort -T ../data/temp_M_rM -o ../data/M_rM.bam



```python
with open('../bin/combos_index.txt', 'w') as fw:
    for combo in compBams:
        command = 'samtools index -b {}'.format(compBams[combo])
        fw.write(command + '\n')
        
```


```bash
%%bash
cat ../bin/combos_index.txt
#parallel --no-notice -j 4 :::: ../bin/combos_index.txt &> ../bin/combos_index.out
```

    samtools index -b ../data/F_rF.bam
    samtools index -b ../data/F_rM.bam
    samtools index -b ../data/M_rF.bam
    samtools index -b ../data/M_rM.bam


#### Generating Coverage Tracks


```python
compBams
```




    {('F', 'rF'): '../data/F_rF.bam',
     ('F', 'rM'): '../data/F_rM.bam',
     ('M', 'rF'): '../data/M_rF.bam',
     ('M', 'rM'): '../data/M_rM.bam'}




```python
bedGraphs = {}
with open('../bin/bedtools_genome_cov.txt', 'w') as fw:
    for combo in compBams:
        bed = compBams[combo].replace('.bam', '.bedgraph')
        command  = 'bedtools genomecov -ibam {} -bga > {}'.format(compBams[combo], bed)
        bedGraphs[combo] = bed
        fw.write(command + '\n')
```


```bash
%%bash
cat ../bin/bedtools_genome_cov.txt
#parallel --no-notice -j 4 :::: ../bin/bedtools_genome_cov.txt &> ../bin/bedtools_genome_cov.out
```

    bedtools genomecov -ibam ../data/F_rF.bam -bga > ../data/F_rF.bedgraph
    bedtools genomecov -ibam ../data/F_rM.bam -bga > ../data/F_rM.bedgraph
    bedtools genomecov -ibam ../data/M_rF.bam -bga > ../data/M_rF.bedgraph
    bedtools genomecov -ibam ../data/M_rM.bam -bga > ../data/M_rM.bedgraph


#### Comparing Coverage Across Male Assembly


```python
m2mCovDF = pd.read_csv('../data/M_rM.bedgraph', sep='\t', names = ['scaffold', 'start', 'stop', 'coverage'])
m2mCovDF = m2mCovDF.merge(maleScaffoldLensDF)
m2mCovDF['interval_length'] = m2mCovDF['stop'] - m2mCovDF['start']
m2mCovDF.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scaffold</th>
      <th>start</th>
      <th>stop</th>
      <th>coverage</th>
      <th>length</th>
      <th>interval_length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>scaffold1</td>
      <td>0</td>
      <td>1</td>
      <td>6</td>
      <td>451</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>scaffold1</td>
      <td>1</td>
      <td>4</td>
      <td>7</td>
      <td>451</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>scaffold1</td>
      <td>4</td>
      <td>5</td>
      <td>8</td>
      <td>451</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>scaffold1</td>
      <td>5</td>
      <td>14</td>
      <td>9</td>
      <td>451</td>
      <td>9</td>
    </tr>
    <tr>
      <th>4</th>
      <td>scaffold1</td>
      <td>14</td>
      <td>19</td>
      <td>10</td>
      <td>451</td>
      <td>5</td>
    </tr>
  </tbody>
</table>
</div>




```python
avgCovList = []
for scaffold in m2mCovDF[m2mCovDF.length > 1000].scaffold.unique():
    df = m2mCovDF[m2mCovDF.scaffold == scaffold]
    avg = sum(df.interval_length * df.coverage) / df.interval_length.sum()
    data = {'scaffold': scaffold, 
            'avg_cov': avg}
    avgCovList.append(data)
    
maleAvgCovDF = pd.DataFrame(avgCovList)
maleAvgCovDF = maleAvgCovDF.merge(maleScaffoldLensDF)
maleAvgCovDF.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov</th>
      <th>scaffold</th>
      <th>length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2361.040028</td>
      <td>scaffold5</td>
      <td>1449</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1797.411190</td>
      <td>scaffold8</td>
      <td>1126</td>
    </tr>
    <tr>
      <th>2</th>
      <td>17.358844</td>
      <td>scaffold9</td>
      <td>1176</td>
    </tr>
    <tr>
      <th>3</th>
      <td>16.534722</td>
      <td>scaffold14</td>
      <td>1152</td>
    </tr>
    <tr>
      <th>4</th>
      <td>27.088109</td>
      <td>scaffold15</td>
      <td>1396</td>
    </tr>
  </tbody>
</table>
</div>




```python
f2mCovDF = pd.read_csv('../data/F_rM.bedgraph', sep='\t', names = ['scaffold', 'start', 'stop', 'coverage'])
f2mCovDF = f2mCovDF.merge(maleScaffoldLensDF)
f2mCovDF['interval_length'] = f2mCovDF['stop'] - f2mCovDF['start']
f2mCovDF.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scaffold</th>
      <th>start</th>
      <th>stop</th>
      <th>coverage</th>
      <th>length</th>
      <th>interval_length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>scaffold1</td>
      <td>0</td>
      <td>1</td>
      <td>7</td>
      <td>451</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>scaffold1</td>
      <td>1</td>
      <td>32</td>
      <td>8</td>
      <td>451</td>
      <td>31</td>
    </tr>
    <tr>
      <th>2</th>
      <td>scaffold1</td>
      <td>32</td>
      <td>36</td>
      <td>9</td>
      <td>451</td>
      <td>4</td>
    </tr>
    <tr>
      <th>3</th>
      <td>scaffold1</td>
      <td>36</td>
      <td>44</td>
      <td>10</td>
      <td>451</td>
      <td>8</td>
    </tr>
    <tr>
      <th>4</th>
      <td>scaffold1</td>
      <td>44</td>
      <td>53</td>
      <td>11</td>
      <td>451</td>
      <td>9</td>
    </tr>
  </tbody>
</table>
</div>




```python
avgCovList = []
for scaffold in f2mCovDF[f2mCovDF.length > 1000].scaffold.unique():
    df = f2mCovDF[f2mCovDF.scaffold == scaffold]
    avg = sum(df.interval_length * df.coverage) / df.interval_length.sum()
    data = {'scaffold': scaffold, 
            'avg_cov': avg}
    avgCovList.append(data)
    
femaleAvgCovDF = pd.DataFrame(avgCovList)
femaleAvgCovDF = femaleAvgCovDF.merge(maleScaffoldLensDF)
femaleAvgCovDF.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov</th>
      <th>scaffold</th>
      <th>length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3155.336784</td>
      <td>scaffold5</td>
      <td>1449</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1574.133215</td>
      <td>scaffold8</td>
      <td>1126</td>
    </tr>
    <tr>
      <th>2</th>
      <td>14.070578</td>
      <td>scaffold9</td>
      <td>1176</td>
    </tr>
    <tr>
      <th>3</th>
      <td>26.012153</td>
      <td>scaffold14</td>
      <td>1152</td>
    </tr>
    <tr>
      <th>4</th>
      <td>24.228510</td>
      <td>scaffold15</td>
      <td>1396</td>
    </tr>
  </tbody>
</table>
</div>




```python
avgCovDf = pd.merge( maleAvgCovDF, femaleAvgCovDF, on='scaffold', how='outer',suffixes=('_male', '_female'))
avgCovDf = avgCovDf.sort_values('avg_cov_male').reset_index(drop=True)
avgCovDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov_male</th>
      <th>scaffold</th>
      <th>length_male</th>
      <th>avg_cov_female</th>
      <th>length_female</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>7.878261</td>
      <td>C6959</td>
      <td>1150</td>
      <td>52.776522</td>
      <td>1150</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11.104535</td>
      <td>C7039</td>
      <td>1301</td>
      <td>42.129900</td>
      <td>1301</td>
    </tr>
    <tr>
      <th>2</th>
      <td>11.339567</td>
      <td>scaffold79</td>
      <td>1016</td>
      <td>19.390748</td>
      <td>1016</td>
    </tr>
    <tr>
      <th>3</th>
      <td>11.996218</td>
      <td>C7049</td>
      <td>1322</td>
      <td>20.648260</td>
      <td>1322</td>
    </tr>
    <tr>
      <th>4</th>
      <td>12.269718</td>
      <td>C6867</td>
      <td>1027</td>
      <td>11.884129</td>
      <td>1027</td>
    </tr>
  </tbody>
</table>
</div>




```python
ax = sns.distplot(avgCovDf.avg_cov_male, kde=False, label='Female')
ax = sns.distplot(avgCovDf.avg_cov_female, kde=False, ax=ax, label='Male')
ax.set_title('Male and Female Average Coverage Distributions')
ax.legend()
```




    <matplotlib.legend.Legend at 0x7f30c6949eb8>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_116_1.png)



```python
avgCovDf['fc'] = np.log2(avgCovDf['avg_cov_male'] / avgCovDf['avg_cov_female'])
```


```python
avgCovDf['fc'].plot()
```




    <matplotlib.axes._subplots.AxesSubplot at 0x7f30c69715f8>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_118_1.png)



```python
avgCovDf.sort_values(['fc'],ascending=False)
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov_male</th>
      <th>scaffold</th>
      <th>length_male</th>
      <th>avg_cov_female</th>
      <th>length_female</th>
      <th>fc</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>77</th>
      <td>36.070813</td>
      <td>C6885</td>
      <td>1045</td>
      <td>0.000000</td>
      <td>1045</td>
      <td>inf</td>
    </tr>
    <tr>
      <th>129</th>
      <td>63.134006</td>
      <td>C7069</td>
      <td>1388</td>
      <td>0.000000</td>
      <td>1388</td>
      <td>inf</td>
    </tr>
    <tr>
      <th>36</th>
      <td>22.005900</td>
      <td>C6859</td>
      <td>1017</td>
      <td>0.000000</td>
      <td>1017</td>
      <td>inf</td>
    </tr>
    <tr>
      <th>59</th>
      <td>30.333333</td>
      <td>C6883</td>
      <td>1044</td>
      <td>0.000000</td>
      <td>1044</td>
      <td>inf</td>
    </tr>
    <tr>
      <th>120</th>
      <td>55.515952</td>
      <td>C7119</td>
      <td>2006</td>
      <td>0.000000</td>
      <td>2006</td>
      <td>inf</td>
    </tr>
    <tr>
      <th>126</th>
      <td>60.450872</td>
      <td>C6925</td>
      <td>1089</td>
      <td>0.071625</td>
      <td>1089</td>
      <td>9.721077</td>
    </tr>
    <tr>
      <th>45</th>
      <td>24.847270</td>
      <td>C6981</td>
      <td>1172</td>
      <td>0.032423</td>
      <td>1172</td>
      <td>9.581845</td>
    </tr>
    <tr>
      <th>32</th>
      <td>21.204601</td>
      <td>C7007</td>
      <td>1217</td>
      <td>0.489729</td>
      <td>1217</td>
      <td>5.436250</td>
    </tr>
    <tr>
      <th>43</th>
      <td>24.285586</td>
      <td>C6941</td>
      <td>1117</td>
      <td>0.841540</td>
      <td>1117</td>
      <td>4.850925</td>
    </tr>
    <tr>
      <th>72</th>
      <td>34.294904</td>
      <td>C6995</td>
      <td>1197</td>
      <td>5.982456</td>
      <td>1197</td>
      <td>2.519184</td>
    </tr>
    <tr>
      <th>6</th>
      <td>13.023969</td>
      <td>C6881</td>
      <td>1043</td>
      <td>2.508150</td>
      <td>1043</td>
      <td>2.376474</td>
    </tr>
    <tr>
      <th>64</th>
      <td>31.868069</td>
      <td>C6889</td>
      <td>1046</td>
      <td>7.540153</td>
      <td>1046</td>
      <td>2.079446</td>
    </tr>
    <tr>
      <th>29</th>
      <td>20.955078</td>
      <td>C6863</td>
      <td>1024</td>
      <td>8.859375</td>
      <td>1024</td>
      <td>1.242023</td>
    </tr>
    <tr>
      <th>52</th>
      <td>27.809568</td>
      <td>C6907</td>
      <td>1066</td>
      <td>12.437148</td>
      <td>1066</td>
      <td>1.160926</td>
    </tr>
    <tr>
      <th>107</th>
      <td>47.389300</td>
      <td>C7005</td>
      <td>1215</td>
      <td>22.535802</td>
      <td>1215</td>
      <td>1.072343</td>
    </tr>
    <tr>
      <th>113</th>
      <td>50.912117</td>
      <td>C7085</td>
      <td>1502</td>
      <td>24.439414</td>
      <td>1502</td>
      <td>1.058799</td>
    </tr>
    <tr>
      <th>39</th>
      <td>23.362117</td>
      <td>C6915</td>
      <td>1077</td>
      <td>11.324977</td>
      <td>1077</td>
      <td>1.044663</td>
    </tr>
    <tr>
      <th>48</th>
      <td>26.098673</td>
      <td>C7003</td>
      <td>1206</td>
      <td>12.917081</td>
      <td>1206</td>
      <td>1.014696</td>
    </tr>
    <tr>
      <th>98</th>
      <td>44.654741</td>
      <td>C6943</td>
      <td>1118</td>
      <td>22.884615</td>
      <td>1118</td>
      <td>0.964435</td>
    </tr>
    <tr>
      <th>42</th>
      <td>24.181818</td>
      <td>C6871</td>
      <td>1034</td>
      <td>12.482592</td>
      <td>1034</td>
      <td>0.954005</td>
    </tr>
    <tr>
      <th>37</th>
      <td>22.266200</td>
      <td>C6955</td>
      <td>1142</td>
      <td>11.534151</td>
      <td>1142</td>
      <td>0.948944</td>
    </tr>
    <tr>
      <th>61</th>
      <td>31.270148</td>
      <td>C7101</td>
      <td>1551</td>
      <td>16.402966</td>
      <td>1551</td>
      <td>0.930829</td>
    </tr>
    <tr>
      <th>104</th>
      <td>46.887690</td>
      <td>C7103</td>
      <td>1576</td>
      <td>28.088198</td>
      <td>1576</td>
      <td>0.739245</td>
    </tr>
    <tr>
      <th>35</th>
      <td>21.975924</td>
      <td>C6971</td>
      <td>1163</td>
      <td>14.024936</td>
      <td>1163</td>
      <td>0.647930</td>
    </tr>
    <tr>
      <th>161</th>
      <td>1949.255493</td>
      <td>scaffold30</td>
      <td>1593</td>
      <td>1257.839297</td>
      <td>1593</td>
      <td>0.631976</td>
    </tr>
    <tr>
      <th>130</th>
      <td>63.540594</td>
      <td>C7125</td>
      <td>2525</td>
      <td>42.431683</td>
      <td>2525</td>
      <td>0.582537</td>
    </tr>
    <tr>
      <th>55</th>
      <td>29.186538</td>
      <td>C6877</td>
      <td>1040</td>
      <td>19.674038</td>
      <td>1040</td>
      <td>0.569010</td>
    </tr>
    <tr>
      <th>86</th>
      <td>39.289428</td>
      <td>C6961</td>
      <td>1154</td>
      <td>26.746967</td>
      <td>1154</td>
      <td>0.554766</td>
    </tr>
    <tr>
      <th>85</th>
      <td>38.777700</td>
      <td>C7079</td>
      <td>1435</td>
      <td>26.468990</td>
      <td>1435</td>
      <td>0.550924</td>
    </tr>
    <tr>
      <th>76</th>
      <td>35.172448</td>
      <td>C6849</td>
      <td>1009</td>
      <td>24.716551</td>
      <td>1009</td>
      <td>0.508968</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>67</th>
      <td>33.095161</td>
      <td>C7015</td>
      <td>1240</td>
      <td>45.112097</td>
      <td>1240</td>
      <td>-0.446894</td>
    </tr>
    <tr>
      <th>53</th>
      <td>28.818526</td>
      <td>C6903</td>
      <td>1058</td>
      <td>39.979206</td>
      <td>1058</td>
      <td>-0.472253</td>
    </tr>
    <tr>
      <th>97</th>
      <td>44.475945</td>
      <td>C6973</td>
      <td>1164</td>
      <td>62.435567</td>
      <td>1164</td>
      <td>-0.489343</td>
    </tr>
    <tr>
      <th>89</th>
      <td>41.233397</td>
      <td>C6899</td>
      <td>1054</td>
      <td>58.067362</td>
      <td>1054</td>
      <td>-0.493914</td>
    </tr>
    <tr>
      <th>40</th>
      <td>23.446729</td>
      <td>C6913</td>
      <td>1070</td>
      <td>33.704673</td>
      <td>1070</td>
      <td>-0.523562</td>
    </tr>
    <tr>
      <th>110</th>
      <td>48.304139</td>
      <td>C6875</td>
      <td>1039</td>
      <td>69.641001</td>
      <td>1039</td>
      <td>-0.527790</td>
    </tr>
    <tr>
      <th>99</th>
      <td>45.140800</td>
      <td>C7017</td>
      <td>1250</td>
      <td>66.401600</td>
      <td>1250</td>
      <td>-0.556786</td>
    </tr>
    <tr>
      <th>23</th>
      <td>18.493056</td>
      <td>C6847</td>
      <td>1008</td>
      <td>27.329365</td>
      <td>1008</td>
      <td>-0.563468</td>
    </tr>
    <tr>
      <th>78</th>
      <td>36.161074</td>
      <td>C6879</td>
      <td>1043</td>
      <td>53.485139</td>
      <td>1043</td>
      <td>-0.564701</td>
    </tr>
    <tr>
      <th>80</th>
      <td>37.267561</td>
      <td>C7021</td>
      <td>1267</td>
      <td>56.493291</td>
      <td>1267</td>
      <td>-0.600159</td>
    </tr>
    <tr>
      <th>7</th>
      <td>13.817606</td>
      <td>C7073</td>
      <td>1420</td>
      <td>21.638028</td>
      <td>1420</td>
      <td>-0.647061</td>
    </tr>
    <tr>
      <th>17</th>
      <td>16.534722</td>
      <td>scaffold14</td>
      <td>1152</td>
      <td>26.012153</td>
      <td>1152</td>
      <td>-0.653687</td>
    </tr>
    <tr>
      <th>15</th>
      <td>16.111024</td>
      <td>C7025</td>
      <td>1270</td>
      <td>25.496850</td>
      <td>1270</td>
      <td>-0.662271</td>
    </tr>
    <tr>
      <th>66</th>
      <td>33.050495</td>
      <td>C6853</td>
      <td>1010</td>
      <td>54.527723</td>
      <td>1010</td>
      <td>-0.722318</td>
    </tr>
    <tr>
      <th>9</th>
      <td>14.180881</td>
      <td>C6911</td>
      <td>1067</td>
      <td>23.728210</td>
      <td>1067</td>
      <td>-0.742656</td>
    </tr>
    <tr>
      <th>14</th>
      <td>15.861406</td>
      <td>C7091</td>
      <td>1508</td>
      <td>26.950928</td>
      <td>1508</td>
      <td>-0.764814</td>
    </tr>
    <tr>
      <th>106</th>
      <td>47.249791</td>
      <td>C6993</td>
      <td>1197</td>
      <td>80.556391</td>
      <td>1197</td>
      <td>-0.769691</td>
    </tr>
    <tr>
      <th>24</th>
      <td>19.157380</td>
      <td>C7121</td>
      <td>2046</td>
      <td>32.688172</td>
      <td>2046</td>
      <td>-0.770868</td>
    </tr>
    <tr>
      <th>2</th>
      <td>11.339567</td>
      <td>scaffold79</td>
      <td>1016</td>
      <td>19.390748</td>
      <td>1016</td>
      <td>-0.774003</td>
    </tr>
    <tr>
      <th>3</th>
      <td>11.996218</td>
      <td>C7049</td>
      <td>1322</td>
      <td>20.648260</td>
      <td>1322</td>
      <td>-0.783441</td>
    </tr>
    <tr>
      <th>57</th>
      <td>29.567111</td>
      <td>C6945</td>
      <td>1125</td>
      <td>51.000000</td>
      <td>1125</td>
      <td>-0.786504</td>
    </tr>
    <tr>
      <th>8</th>
      <td>14.152971</td>
      <td>C7127</td>
      <td>2726</td>
      <td>25.190389</td>
      <td>2726</td>
      <td>-0.831768</td>
    </tr>
    <tr>
      <th>34</th>
      <td>21.775021</td>
      <td>C6977</td>
      <td>1169</td>
      <td>40.591959</td>
      <td>1169</td>
      <td>-0.898520</td>
    </tr>
    <tr>
      <th>33</th>
      <td>21.741935</td>
      <td>C7013</td>
      <td>1240</td>
      <td>40.694355</td>
      <td>1240</td>
      <td>-0.904348</td>
    </tr>
    <tr>
      <th>11</th>
      <td>14.848562</td>
      <td>scaffold20</td>
      <td>1182</td>
      <td>28.628596</td>
      <td>1182</td>
      <td>-0.947134</td>
    </tr>
    <tr>
      <th>25</th>
      <td>19.452290</td>
      <td>C6893</td>
      <td>1048</td>
      <td>41.033397</td>
      <td>1048</td>
      <td>-1.076859</td>
    </tr>
    <tr>
      <th>5</th>
      <td>12.605697</td>
      <td>C7059</td>
      <td>1334</td>
      <td>29.838831</td>
      <td>1334</td>
      <td>-1.243115</td>
    </tr>
    <tr>
      <th>22</th>
      <td>17.769287</td>
      <td>C7065</td>
      <td>1361</td>
      <td>45.359295</td>
      <td>1361</td>
      <td>-1.352012</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11.104535</td>
      <td>C7039</td>
      <td>1301</td>
      <td>42.129900</td>
      <td>1301</td>
      <td>-1.923696</td>
    </tr>
    <tr>
      <th>0</th>
      <td>7.878261</td>
      <td>C6959</td>
      <td>1150</td>
      <td>52.776522</td>
      <td>1150</td>
      <td>-2.743947</td>
    </tr>
  </tbody>
</table>
<p>168 rows × 6 columns</p>
</div>



### Aligning Individual Libararies to Male Genome

#### Converting to FQ


```python
unmappedFastqs={}
with open('../bin/extract_indvidual_fq.txt','w') as fw:
    for sample_name in unmappedPathsDict:
        sex = sampleToSexDict[sample_name]
        if sex =='M' or sex =='F':
            key = (sample_name, sex)
            output_1 = '../data/unmapped_{}_{}_1.fq'.format(*key)
            output_2 = '../data/unmapped_{}_{}_2.fq'.format(*key)
            command = "bedtools bamtofastq -i {} -fq {} -fq2 {}".format(
                unmappedPathsDict[sample_name],
                output_1,
                output_2
                )
            fw.write(command + '\n')
            unmappedFastqs[key] = [output_1, output_2]
```


```bash
%%bash
cat ../bin/extract_indvidual_fq.txt
#parallel --no-notice -j 8 :::: ../bin/extract_indvidual_fq.txt &> ../bin/extract_indvidual_fq.out
```

    bedtools bamtofastq -i ../data/ATCACG_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_21545_M_1.fq -fq2 ../data/unmapped_A.tig_21545_M_2.fq
    bedtools bamtofastq -i ../data/ACAGTG_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_21321_F_1.fq -fq2 ../data/unmapped_A.tig_21321_F_2.fq
    bedtools bamtofastq -i ../data/GCCAAT_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_21322_M_1.fq -fq2 ../data/unmapped_A.tig_21322_M_2.fq
    bedtools bamtofastq -i ../data/CTTGTA_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_21544_F_1.fq -fq2 ../data/unmapped_A.tig_21544_F_2.fq
    bedtools bamtofastq -i ../data/CGATGT_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_17115_M_1.fq -fq2 ../data/unmapped_A.tig_17115_M_2.fq
    bedtools bamtofastq -i ../data/TTAGGC_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_21546_F_1.fq -fq2 ../data/unmapped_A.tig_21546_F_2.fq
    bedtools bamtofastq -i ../data/TGACCA_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_17118_F_1.fq -fq2 ../data/unmapped_A.tig_17118_F_2.fq
    bedtools bamtofastq -i ../data/CAGATC_merged_sorted_markDup_indelRealign.unmapped.bam -fq ../data/unmapped_A.tig_21323_M_1.fq -fq2 ../data/unmapped_A.tig_21323_M_2.fq


#### Aligning


```python
indvBams = {}
with open('../bin/indv_to_male_align.txt', 'w') as fw:
    for key in unmappedFastqs:
        out_bam = '../data/{}_{}_rM.bam'.format(*key)
        temp = '../data/temp_{}_{}_rM_'.format(*key)
        stdout = '../data/{}_{}_rM.out'.format(*key)
        bwa_command = 'bwa mem -t 8 {} {} {} | samtools sort -T {} -o {} &> {}'.format(
            assemblyPaths['rM'],
            unmappedFastqs[key][0],
            unmappedFastqs[key][1],
            temp,
            out_bam,
            stdout,
        )
        indvBams[key] = out_bam
        fw.write(bwa_command + '\n')
```


```bash
%%bash
cat ../bin/indv_to_male_align.txt
#parallel --no-notice -j 4 :::: ../bin/indv_to_male_align.txt &> ../bin/indv_to_male_align.out
```

    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_21545_M_1.fq ../data/unmapped_A.tig_21545_M_2.fq | samtools sort -T ../data/temp_A.tig_21545_M_rM_ -o ../data/A.tig_21545_M_rM.bam &> ../data/A.tig_21545_M_rM.out
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_21321_F_1.fq ../data/unmapped_A.tig_21321_F_2.fq | samtools sort -T ../data/temp_A.tig_21321_F_rM_ -o ../data/A.tig_21321_F_rM.bam &> ../data/A.tig_21321_F_rM.out
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_21322_M_1.fq ../data/unmapped_A.tig_21322_M_2.fq | samtools sort -T ../data/temp_A.tig_21322_M_rM_ -o ../data/A.tig_21322_M_rM.bam &> ../data/A.tig_21322_M_rM.out
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_21544_F_1.fq ../data/unmapped_A.tig_21544_F_2.fq | samtools sort -T ../data/temp_A.tig_21544_F_rM_ -o ../data/A.tig_21544_F_rM.bam &> ../data/A.tig_21544_F_rM.out
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_17115_M_1.fq ../data/unmapped_A.tig_17115_M_2.fq | samtools sort -T ../data/temp_A.tig_17115_M_rM_ -o ../data/A.tig_17115_M_rM.bam &> ../data/A.tig_17115_M_rM.out
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_21546_F_1.fq ../data/unmapped_A.tig_21546_F_2.fq | samtools sort -T ../data/temp_A.tig_21546_F_rM_ -o ../data/A.tig_21546_F_rM.bam &> ../data/A.tig_21546_F_rM.out
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_17118_F_1.fq ../data/unmapped_A.tig_17118_F_2.fq | samtools sort -T ../data/temp_A.tig_17118_F_rM_ -o ../data/A.tig_17118_F_rM.bam &> ../data/A.tig_17118_F_rM.out
    bwa mem -t 8 ../data/male_assembly/male_assembly.scafSeq ../data/unmapped_A.tig_21323_M_1.fq ../data/unmapped_A.tig_21323_M_2.fq | samtools sort -T ../data/temp_A.tig_21323_M_rM_ -o ../data/A.tig_21323_M_rM.bam &> ../data/A.tig_21323_M_rM.out



```python
with open('../bin/indv_to_male_index.txt', 'w') as fw:
    for key in indvBams:
        command = 'samtools index -b {}'.format(indvBams[key])
        fw.write(command + '\n')
        

```


```bash
%%bash
cat ../bin/indv_to_male_index.txt
#parallel --no-notice -j 8 :::: ../bin/indv_to_male_index.txt &> ../bin/indv_to_male_index.out
```

    samtools index -b ../data/A.tig_21545_M_rM.bam
    samtools index -b ../data/A.tig_21321_F_rM.bam
    samtools index -b ../data/A.tig_21322_M_rM.bam
    samtools index -b ../data/A.tig_21544_F_rM.bam
    samtools index -b ../data/A.tig_17115_M_rM.bam
    samtools index -b ../data/A.tig_21546_F_rM.bam
    samtools index -b ../data/A.tig_17118_F_rM.bam
    samtools index -b ../data/A.tig_21323_M_rM.bam



```bash
%%bash
rm ../data/flagstat_indv_data.csv
touch ../data/flagstat_indv_data.csv
for file in ../data/A.tig_*_*_rM.bam
 do
  sample=$(basename $file | sed 's/.bam//')
  flagstat=$(echo $file | sed 's/.bam/.flagstat/')
  echo $file
  echo $sample
  samtools flagstat $file > $flagstat 
  grep "mapped (" $flagstat | sed "s/ + 0 mapped (/,/" | sed "s/%.*$//" | sed "s/^/$sample,/" >> ../data/flagstat_indv_data.csv
 done

head -n 5 ../data/flagstat_indv_data.csv
```

    ../data/A.tig_17115_M_rM.bam
    A.tig_17115_M_rM
    ../data/A.tig_17118_F_rM.bam
    A.tig_17118_F_rM
    ../data/A.tig_21321_F_rM.bam
    A.tig_21321_F_rM
    ../data/A.tig_21322_M_rM.bam
    A.tig_21322_M_rM
    ../data/A.tig_21323_M_rM.bam
    A.tig_21323_M_rM
    ../data/A.tig_21544_F_rM.bam
    A.tig_21544_F_rM
    ../data/A.tig_21545_M_rM.bam
    A.tig_21545_M_rM
    ../data/A.tig_21546_F_rM.bam
    A.tig_21546_F_rM
    A.tig_17115_M_rM,724129,10.20
    A.tig_17118_F_rM,543842,11.00
    A.tig_21321_F_rM,749486,12.00
    A.tig_21322_M_rM,774808,12.94
    A.tig_21323_M_rM,793697,14.20



```python
normFactorsDF = pd.read_csv('../data/flagstat_indv_data.csv', names = ['alignment','mapped_reads','percent_mapped'])
normFactorsDF['sample_name'] = normFactorsDF['alignment'].apply(lambda x: '_'.join(x.split('_')[:2]))
normFactorsDF['sex'] = normFactorsDF['alignment'].apply(lambda x: x.split('_')[2])


normFactorsDF
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>alignment</th>
      <th>mapped_reads</th>
      <th>percent_mapped</th>
      <th>sample_name</th>
      <th>sex</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>A.tig_17115_M_rM</td>
      <td>724129</td>
      <td>10.20</td>
      <td>A.tig_17115</td>
      <td>M</td>
    </tr>
    <tr>
      <th>1</th>
      <td>A.tig_17118_F_rM</td>
      <td>543842</td>
      <td>11.00</td>
      <td>A.tig_17118</td>
      <td>F</td>
    </tr>
    <tr>
      <th>2</th>
      <td>A.tig_21321_F_rM</td>
      <td>749486</td>
      <td>12.00</td>
      <td>A.tig_21321</td>
      <td>F</td>
    </tr>
    <tr>
      <th>3</th>
      <td>A.tig_21322_M_rM</td>
      <td>774808</td>
      <td>12.94</td>
      <td>A.tig_21322</td>
      <td>M</td>
    </tr>
    <tr>
      <th>4</th>
      <td>A.tig_21323_M_rM</td>
      <td>793697</td>
      <td>14.20</td>
      <td>A.tig_21323</td>
      <td>M</td>
    </tr>
    <tr>
      <th>5</th>
      <td>A.tig_21544_F_rM</td>
      <td>964352</td>
      <td>14.41</td>
      <td>A.tig_21544</td>
      <td>F</td>
    </tr>
    <tr>
      <th>6</th>
      <td>A.tig_21545_M_rM</td>
      <td>851393</td>
      <td>13.71</td>
      <td>A.tig_21545</td>
      <td>M</td>
    </tr>
    <tr>
      <th>7</th>
      <td>A.tig_21546_F_rM</td>
      <td>1271911</td>
      <td>18.61</td>
      <td>A.tig_21546</td>
      <td>F</td>
    </tr>
  </tbody>
</table>
</div>



#### Generating Coverage Tracks


```python
indvBedGraphs = {}
with open('../bin/indv_bedtools_genome_cov.txt', 'w') as fw:
    for key in indvBams:
        bed = indvBams[key].replace('.bam', '.bedgraph')
        command  = 'bedtools genomecov -ibam {} -bga > {}'.format(indvBams[key], bed)
        indvBedGraphs[key] = bed
        fw.write(command + '\n')
```


```bash
%%bash
cat ../bin/indv_bedtools_genome_cov.txt
#parallel --no-notice -j 8 :::: ../bin/indv_bedtools_genome_cov.txt &> ../bin/indv_bedtools_genome_cov.out
```

    bedtools genomecov -ibam ../data/A.tig_21545_M_rM.bam -bga > ../data/A.tig_21545_M_rM.bedgraph
    bedtools genomecov -ibam ../data/A.tig_21321_F_rM.bam -bga > ../data/A.tig_21321_F_rM.bedgraph
    bedtools genomecov -ibam ../data/A.tig_21322_M_rM.bam -bga > ../data/A.tig_21322_M_rM.bedgraph
    bedtools genomecov -ibam ../data/A.tig_21544_F_rM.bam -bga > ../data/A.tig_21544_F_rM.bedgraph
    bedtools genomecov -ibam ../data/A.tig_17115_M_rM.bam -bga > ../data/A.tig_17115_M_rM.bedgraph
    bedtools genomecov -ibam ../data/A.tig_21546_F_rM.bam -bga > ../data/A.tig_21546_F_rM.bedgraph
    bedtools genomecov -ibam ../data/A.tig_17118_F_rM.bam -bga > ../data/A.tig_17118_F_rM.bedgraph
    bedtools genomecov -ibam ../data/A.tig_21323_M_rM.bam -bga > ../data/A.tig_21323_M_rM.bedgraph


#### Load In Coverage Files


```python
allSampList = []
for key in indvBedGraphs:
    mapped_reads = normFactorsDF[normFactorsDF['sample_name']==key[0]]['mapped_reads'].iloc[0]
    df = pd.read_csv(indvBedGraphs[key], sep='\t', names = ['scaffold', 'start', 'stop', 'coverage'])
    df['interval_length'] = df.stop - df.start
    for scaffold in df.scaffold.unique():
        sdf = df[df.scaffold == scaffold]
        
        avg = ((sum(sdf.interval_length * sdf.coverage) / sdf.interval_length.sum()) / mapped_reads) * 1e6
        
        data = {'scaffold': scaffold, 'avg_cov': avg, 'sample_name':key[0], 'sex':key[1]}
        
        allSampList.append(data)

allSampsAvgDf = pd.DataFrame(allSampList)
allSampsAvgDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov</th>
      <th>sample_name</th>
      <th>scaffold</th>
      <th>sex</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3.586141</td>
      <td>A.tig_21545</td>
      <td>scaffold1</td>
      <td>M</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3.778757</td>
      <td>A.tig_21545</td>
      <td>scaffold2</td>
      <td>M</td>
    </tr>
    <tr>
      <th>2</th>
      <td>4.216106</td>
      <td>A.tig_21545</td>
      <td>scaffold3</td>
      <td>M</td>
    </tr>
    <tr>
      <th>3</th>
      <td>6.164294</td>
      <td>A.tig_21545</td>
      <td>scaffold4</td>
      <td>M</td>
    </tr>
    <tr>
      <th>4</th>
      <td>862.024958</td>
      <td>A.tig_21545</td>
      <td>scaffold5</td>
      <td>M</td>
    </tr>
  </tbody>
</table>
</div>




```python
allSampsAvgDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov</th>
      <th>sample_name</th>
      <th>scaffold</th>
      <th>sex</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3.586141</td>
      <td>A.tig_21545</td>
      <td>scaffold1</td>
      <td>M</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3.778757</td>
      <td>A.tig_21545</td>
      <td>scaffold2</td>
      <td>M</td>
    </tr>
    <tr>
      <th>2</th>
      <td>4.216106</td>
      <td>A.tig_21545</td>
      <td>scaffold3</td>
      <td>M</td>
    </tr>
    <tr>
      <th>3</th>
      <td>6.164294</td>
      <td>A.tig_21545</td>
      <td>scaffold4</td>
      <td>M</td>
    </tr>
    <tr>
      <th>4</th>
      <td>862.024958</td>
      <td>A.tig_21545</td>
      <td>scaffold5</td>
      <td>M</td>
    </tr>
  </tbody>
</table>
</div>




```python
allSampsAvgDf = allSampsAvgDf.merge(maleScaffoldLensDF)
allSampsAvgDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov</th>
      <th>sample_name</th>
      <th>scaffold</th>
      <th>sex</th>
      <th>length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3.586141</td>
      <td>A.tig_21545</td>
      <td>scaffold1</td>
      <td>M</td>
      <td>451</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.106503</td>
      <td>A.tig_21321</td>
      <td>scaffold1</td>
      <td>F</td>
      <td>451</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3.316750</td>
      <td>A.tig_21322</td>
      <td>scaffold1</td>
      <td>M</td>
      <td>451</td>
    </tr>
    <tr>
      <th>3</th>
      <td>2.805096</td>
      <td>A.tig_21544</td>
      <td>scaffold1</td>
      <td>F</td>
      <td>451</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1.674923</td>
      <td>A.tig_17115</td>
      <td>scaffold1</td>
      <td>M</td>
      <td>451</td>
    </tr>
  </tbody>
</table>
</div>




```python
allSampsAvgDf = allSampsAvgDf.sort_values('length', ascending=False)
allSampsAvgDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov</th>
      <th>sample_name</th>
      <th>scaffold</th>
      <th>sex</th>
      <th>length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>267</th>
      <td>1226.652593</td>
      <td>A.tig_21544</td>
      <td>scaffold34</td>
      <td>F</td>
      <td>4646</td>
    </tr>
    <tr>
      <th>265</th>
      <td>1073.736011</td>
      <td>A.tig_21321</td>
      <td>scaffold34</td>
      <td>F</td>
      <td>4646</td>
    </tr>
    <tr>
      <th>266</th>
      <td>1076.089569</td>
      <td>A.tig_21322</td>
      <td>scaffold34</td>
      <td>M</td>
      <td>4646</td>
    </tr>
    <tr>
      <th>268</th>
      <td>760.383873</td>
      <td>A.tig_17115</td>
      <td>scaffold34</td>
      <td>M</td>
      <td>4646</td>
    </tr>
    <tr>
      <th>269</th>
      <td>1247.072896</td>
      <td>A.tig_21546</td>
      <td>scaffold34</td>
      <td>F</td>
      <td>4646</td>
    </tr>
  </tbody>
</table>
</div>




```python
allSampsAvgDf['sex_name'] = allSampsAvgDf.apply(lambda x: '_'.join([x['sex'], x['sample_name']]), axis=1)
```

#### Clustering and Identification of >=1Kb Shared Male Contigs/Scaffolds


```python
pivAvgCov = allSampsAvgDf[allSampsAvgDf['length'] > 1000].pivot(index='scaffold', columns='sex_name', values='avg_cov').reset_index()
pivAvgCov = pivAvgCov.merge(maleScaffoldLensDF)
pivAvgCov.sort_values('length', ascending=False, inplace=True)
pivAvgCov = pivAvgCov.set_index('scaffold')
pivAvgCov
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>F_A.tig_17118</th>
      <th>F_A.tig_21321</th>
      <th>F_A.tig_21544</th>
      <th>F_A.tig_21546</th>
      <th>M_A.tig_17115</th>
      <th>M_A.tig_21322</th>
      <th>M_A.tig_21323</th>
      <th>M_A.tig_21545</th>
      <th>length</th>
    </tr>
    <tr>
      <th>scaffold</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>scaffold34</th>
      <td>1395.921579</td>
      <td>1073.736011</td>
      <td>1226.652593</td>
      <td>1247.072896</td>
      <td>760.383873</td>
      <td>1076.089569</td>
      <td>1286.891138</td>
      <td>1147.491726</td>
      <td>4646</td>
    </tr>
    <tr>
      <th>scaffold36</th>
      <td>810.090376</td>
      <td>631.046307</td>
      <td>731.777908</td>
      <td>695.487032</td>
      <td>444.348945</td>
      <td>616.281025</td>
      <td>755.770015</td>
      <td>666.425368</td>
      <td>3760</td>
    </tr>
    <tr>
      <th>C7127</th>
      <td>0.000000</td>
      <td>10.618674</td>
      <td>0.000000</td>
      <td>13.547993</td>
      <td>19.544821</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>2726</td>
    </tr>
    <tr>
      <th>C7125</th>
      <td>17.339049</td>
      <td>12.366495</td>
      <td>24.610791</td>
      <td>0.000000</td>
      <td>20.922917</td>
      <td>15.052212</td>
      <td>29.413923</td>
      <td>15.717050</td>
      <td>2525</td>
    </tr>
    <tr>
      <th>C7123</th>
      <td>37.711805</td>
      <td>14.268528</td>
      <td>19.764912</td>
      <td>10.023259</td>
      <td>19.267968</td>
      <td>19.728127</td>
      <td>15.996579</td>
      <td>9.289792</td>
      <td>2105</td>
    </tr>
    <tr>
      <th>C7121</th>
      <td>18.235811</td>
      <td>8.888464</td>
      <td>0.000000</td>
      <td>12.665189</td>
      <td>0.000000</td>
      <td>0.316037</td>
      <td>14.907285</td>
      <td>8.316541</td>
      <td>2046</td>
    </tr>
    <tr>
      <th>C7119</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.461638</td>
      <td>64.775669</td>
      <td>2006</td>
    </tr>
    <tr>
      <th>C7117</th>
      <td>38.956946</td>
      <td>21.788658</td>
      <td>25.149309</td>
      <td>26.448646</td>
      <td>37.753493</td>
      <td>24.135142</td>
      <td>27.286284</td>
      <td>26.797969</td>
      <td>1974</td>
    </tr>
    <tr>
      <th>C7115</th>
      <td>57.053511</td>
      <td>27.361103</td>
      <td>31.870241</td>
      <td>36.641792</td>
      <td>51.255377</td>
      <td>35.665242</td>
      <td>36.449376</td>
      <td>41.716539</td>
      <td>1922</td>
    </tr>
    <tr>
      <th>scaffold33</th>
      <td>1361.163440</td>
      <td>990.958659</td>
      <td>1175.303032</td>
      <td>1237.315382</td>
      <td>726.564855</td>
      <td>998.046332</td>
      <td>1117.411570</td>
      <td>1112.341537</td>
      <td>1784</td>
    </tr>
    <tr>
      <th>C7113</th>
      <td>33.188159</td>
      <td>13.403437</td>
      <td>22.281439</td>
      <td>24.943019</td>
      <td>28.374150</td>
      <td>18.352772</td>
      <td>29.769814</td>
      <td>22.058116</td>
      <td>1751</td>
    </tr>
    <tr>
      <th>C7111</th>
      <td>18.989512</td>
      <td>18.445527</td>
      <td>22.224764</td>
      <td>23.355883</td>
      <td>29.052820</td>
      <td>1.089239</td>
      <td>24.963631</td>
      <td>21.406113</td>
      <td>1711</td>
    </tr>
    <tr>
      <th>C7109</th>
      <td>30.044256</td>
      <td>14.449551</td>
      <td>27.176212</td>
      <td>24.505225</td>
      <td>33.077032</td>
      <td>19.423858</td>
      <td>22.473979</td>
      <td>24.159870</td>
      <td>1668</td>
    </tr>
    <tr>
      <th>C7107</th>
      <td>0.000000</td>
      <td>8.975922</td>
      <td>12.266164</td>
      <td>0.214855</td>
      <td>15.171460</td>
      <td>7.341991</td>
      <td>0.000000</td>
      <td>10.917452</td>
      <td>1654</td>
    </tr>
    <tr>
      <th>C7105</th>
      <td>11.554877</td>
      <td>7.069274</td>
      <td>9.652998</td>
      <td>12.234518</td>
      <td>12.181541</td>
      <td>7.570592</td>
      <td>20.225386</td>
      <td>8.920416</td>
      <td>1609</td>
    </tr>
    <tr>
      <th>scaffold30</th>
      <td>803.864250</td>
      <td>240.279340</td>
      <td>271.119351</td>
      <td>296.616118</td>
      <td>430.152467</td>
      <td>1008.178002</td>
      <td>569.434644</td>
      <td>466.731878</td>
      <td>1593</td>
    </tr>
    <tr>
      <th>C7103</th>
      <td>0.000000</td>
      <td>6.611130</td>
      <td>11.684946</td>
      <td>9.328363</td>
      <td>32.470306</td>
      <td>8.350685</td>
      <td>11.416886</td>
      <td>9.212284</td>
      <td>1576</td>
    </tr>
    <tr>
      <th>C7101</th>
      <td>0.000000</td>
      <td>7.873868</td>
      <td>10.889812</td>
      <td>0.000000</td>
      <td>17.799458</td>
      <td>0.000000</td>
      <td>13.209328</td>
      <td>9.275200</td>
      <td>1551</td>
    </tr>
    <tr>
      <th>C7099</th>
      <td>30.677476</td>
      <td>6.524523</td>
      <td>9.437327</td>
      <td>7.897275</td>
      <td>19.159387</td>
      <td>16.018658</td>
      <td>16.015575</td>
      <td>11.832067</td>
      <td>1546</td>
    </tr>
    <tr>
      <th>C7097</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>8.553273</td>
      <td>11.007573</td>
      <td>18.286562</td>
      <td>0.000000</td>
      <td>0.238810</td>
      <td>12.447114</td>
      <td>1530</td>
    </tr>
    <tr>
      <th>C7095</th>
      <td>0.000000</td>
      <td>6.100672</td>
      <td>11.592868</td>
      <td>9.129445</td>
      <td>12.062586</td>
      <td>0.000000</td>
      <td>17.452471</td>
      <td>0.000000</td>
      <td>1520</td>
    </tr>
    <tr>
      <th>C7093</th>
      <td>30.272719</td>
      <td>13.339827</td>
      <td>21.599791</td>
      <td>22.883645</td>
      <td>27.947711</td>
      <td>20.198125</td>
      <td>18.918091</td>
      <td>21.776544</td>
      <td>1510</td>
    </tr>
    <tr>
      <th>scaffold52</th>
      <td>1127.599409</td>
      <td>558.233692</td>
      <td>528.212077</td>
      <td>536.542828</td>
      <td>718.081216</td>
      <td>421.090662</td>
      <td>817.505367</td>
      <td>658.493605</td>
      <td>1509</td>
    </tr>
    <tr>
      <th>C7091</th>
      <td>15.051571</td>
      <td>7.065850</td>
      <td>0.000000</td>
      <td>10.589967</td>
      <td>0.000000</td>
      <td>9.585672</td>
      <td>0.000000</td>
      <td>9.906530</td>
      <td>1508</td>
    </tr>
    <tr>
      <th>C7089</th>
      <td>24.356058</td>
      <td>17.364725</td>
      <td>17.307337</td>
      <td>19.021264</td>
      <td>21.572486</td>
      <td>21.109077</td>
      <td>24.294400</td>
      <td>17.101542</td>
      <td>1505</td>
    </tr>
    <tr>
      <th>C7085</th>
      <td>3.662848</td>
      <td>3.128642</td>
      <td>2.435696</td>
      <td>14.176012</td>
      <td>24.457502</td>
      <td>13.566352</td>
      <td>22.779340</td>
      <td>5.197873</td>
      <td>1502</td>
    </tr>
    <tr>
      <th>C7087</th>
      <td>29.891632</td>
      <td>2.852376</td>
      <td>10.762489</td>
      <td>8.366266</td>
      <td>14.674869</td>
      <td>18.901637</td>
      <td>10.988708</td>
      <td>11.645363</td>
      <td>1502</td>
    </tr>
    <tr>
      <th>scaffold61</th>
      <td>17.843568</td>
      <td>10.106700</td>
      <td>11.085658</td>
      <td>6.194011</td>
      <td>9.909160</td>
      <td>10.945876</td>
      <td>4.765437</td>
      <td>11.782212</td>
      <td>1470</td>
    </tr>
    <tr>
      <th>scaffold5</th>
      <td>1025.010489</td>
      <td>780.622403</td>
      <td>914.961617</td>
      <td>887.074255</td>
      <td>502.367890</td>
      <td>749.791044</td>
      <td>856.763163</td>
      <td>862.024958</td>
      <td>1449</td>
    </tr>
    <tr>
      <th>C7083</th>
      <td>0.000000</td>
      <td>14.708097</td>
      <td>11.751561</td>
      <td>9.227465</td>
      <td>28.804577</td>
      <td>7.146731</td>
      <td>12.660259</td>
      <td>11.102950</td>
      <td>1446</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>C6893</th>
      <td>28.211424</td>
      <td>11.257079</td>
      <td>5.336218</td>
      <td>9.519396</td>
      <td>0.000000</td>
      <td>7.269715</td>
      <td>6.214276</td>
      <td>10.438663</td>
      <td>1048</td>
    </tr>
    <tr>
      <th>C6891</th>
      <td>35.465242</td>
      <td>13.059571</td>
      <td>18.931806</td>
      <td>24.817322</td>
      <td>35.462028</td>
      <td>17.879156</td>
      <td>20.786985</td>
      <td>21.210254</td>
      <td>1047</td>
    </tr>
    <tr>
      <th>C6889</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>7.818880</td>
      <td>0.000000</td>
      <td>13.564130</td>
      <td>7.161461</td>
      <td>12.357158</td>
      <td>7.856880</td>
      <td>1046</td>
    </tr>
    <tr>
      <th>C6887</th>
      <td>28.820290</td>
      <td>22.899266</td>
      <td>38.255602</td>
      <td>26.217567</td>
      <td>43.883111</td>
      <td>18.020825</td>
      <td>34.919862</td>
      <td>38.717298</td>
      <td>1045</td>
    </tr>
    <tr>
      <th>C6885</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>13.512356</td>
      <td>10.490637</td>
      <td>10.605086</td>
      <td>11.733094</td>
      <td>1045</td>
    </tr>
    <tr>
      <th>C6883</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>13.303074</td>
      <td>6.477937</td>
      <td>11.219864</td>
      <td>7.666048</td>
      <td>1044</td>
    </tr>
    <tr>
      <th>C6879</th>
      <td>51.173497</td>
      <td>5.416304</td>
      <td>10.906533</td>
      <td>8.674033</td>
      <td>6.608263</td>
      <td>14.034962</td>
      <td>11.100159</td>
      <td>13.686892</td>
      <td>1043</td>
    </tr>
    <tr>
      <th>C6881</th>
      <td>1.018992</td>
      <td>0.000000</td>
      <td>0.960411</td>
      <td>0.593998</td>
      <td>17.299853</td>
      <td>0.617479</td>
      <td>0.788813</td>
      <td>0.042793</td>
      <td>1043</td>
    </tr>
    <tr>
      <th>C6877</th>
      <td>10.376670</td>
      <td>4.177222</td>
      <td>6.899811</td>
      <td>3.338405</td>
      <td>17.430755</td>
      <td>12.654500</td>
      <td>5.066359</td>
      <td>3.216448</td>
      <td>1040</td>
    </tr>
    <tr>
      <th>C6873</th>
      <td>9.291183</td>
      <td>0.056503</td>
      <td>12.338795</td>
      <td>0.000000</td>
      <td>9.533873</td>
      <td>4.973755</td>
      <td>5.529611</td>
      <td>6.875445</td>
      <td>1039</td>
    </tr>
    <tr>
      <th>C6875</th>
      <td>34.380916</td>
      <td>17.150027</td>
      <td>8.180951</td>
      <td>23.743951</td>
      <td>13.351143</td>
      <td>24.517235</td>
      <td>10.661477</td>
      <td>13.129138</td>
      <td>1039</td>
    </tr>
    <tr>
      <th>C6871</th>
      <td>0.000000</td>
      <td>8.079038</td>
      <td>6.624948</td>
      <td>0.000000</td>
      <td>14.035404</td>
      <td>0.000000</td>
      <td>9.132640</td>
      <td>7.951470</td>
      <td>1034</td>
    </tr>
    <tr>
      <th>C6869</th>
      <td>11.214531</td>
      <td>7.472305</td>
      <td>9.022910</td>
      <td>5.578262</td>
      <td>15.513470</td>
      <td>6.145260</td>
      <td>6.023451</td>
      <td>10.447875</td>
      <td>1031</td>
    </tr>
    <tr>
      <th>C6867</th>
      <td>0.000000</td>
      <td>7.467630</td>
      <td>6.519657</td>
      <td>0.000000</td>
      <td>10.734449</td>
      <td>5.803492</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1027</td>
    </tr>
    <tr>
      <th>C6865</th>
      <td>41.741498</td>
      <td>20.609316</td>
      <td>30.840130</td>
      <td>31.689358</td>
      <td>31.828249</td>
      <td>22.832796</td>
      <td>28.961121</td>
      <td>32.641152</td>
      <td>1026</td>
    </tr>
    <tr>
      <th>scaffold76</th>
      <td>6.077258</td>
      <td>8.021092</td>
      <td>7.897515</td>
      <td>11.284458</td>
      <td>6.974838</td>
      <td>6.683414</td>
      <td>12.489975</td>
      <td>8.822830</td>
      <td>1026</td>
    </tr>
    <tr>
      <th>C6863</th>
      <td>6.214825</td>
      <td>0.000000</td>
      <td>2.764567</td>
      <td>2.212007</td>
      <td>7.612864</td>
      <td>0.620113</td>
      <td>9.770584</td>
      <td>8.464988</td>
      <td>1024</td>
    </tr>
    <tr>
      <th>C6861</th>
      <td>41.232248</td>
      <td>18.341997</td>
      <td>28.178329</td>
      <td>27.071562</td>
      <td>41.414248</td>
      <td>25.919980</td>
      <td>31.562147</td>
      <td>30.046118</td>
      <td>1024</td>
    </tr>
    <tr>
      <th>C6859</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>12.770912</td>
      <td>7.035714</td>
      <td>4.557788</td>
      <td>4.333231</td>
      <td>1017</td>
    </tr>
    <tr>
      <th>scaffold79</th>
      <td>9.298816</td>
      <td>3.741409</td>
      <td>5.350172</td>
      <td>5.008274</td>
      <td>3.517666</td>
      <td>3.044950</td>
      <td>3.957112</td>
      <td>3.866984</td>
      <td>1016</td>
    </tr>
    <tr>
      <th>C6857</th>
      <td>0.000000</td>
      <td>5.469758</td>
      <td>0.000000</td>
      <td>9.203792</td>
      <td>11.107620</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>10.362618</td>
      <td>1015</td>
    </tr>
    <tr>
      <th>scaffold69</th>
      <td>12.203402</td>
      <td>11.802758</td>
      <td>11.931760</td>
      <td>13.384342</td>
      <td>18.405200</td>
      <td>15.314433</td>
      <td>13.887800</td>
      <td>12.546653</td>
      <td>1013</td>
    </tr>
    <tr>
      <th>C6855</th>
      <td>25.069342</td>
      <td>17.616548</td>
      <td>23.949507</td>
      <td>22.608634</td>
      <td>32.550253</td>
      <td>20.550899</td>
      <td>24.524416</td>
      <td>25.111857</td>
      <td>1013</td>
    </tr>
    <tr>
      <th>C6853</th>
      <td>14.227705</td>
      <td>11.814037</td>
      <td>14.412797</td>
      <td>18.898047</td>
      <td>18.677269</td>
      <td>6.646169</td>
      <td>7.907599</td>
      <td>9.513820</td>
      <td>1010</td>
    </tr>
    <tr>
      <th>C6851</th>
      <td>20.521394</td>
      <td>10.230113</td>
      <td>13.146878</td>
      <td>18.656732</td>
      <td>31.838865</td>
      <td>16.484442</td>
      <td>21.501085</td>
      <td>20.760386</td>
      <td>1010</td>
    </tr>
    <tr>
      <th>C6849</th>
      <td>0.072895</td>
      <td>3.718439</td>
      <td>8.452967</td>
      <td>10.801349</td>
      <td>24.802703</td>
      <td>6.382859</td>
      <td>8.041554</td>
      <td>6.911078</td>
      <td>1009</td>
    </tr>
    <tr>
      <th>scaffold70</th>
      <td>23.210815</td>
      <td>12.590639</td>
      <td>11.221451</td>
      <td>18.390338</td>
      <td>16.462032</td>
      <td>13.597839</td>
      <td>14.562902</td>
      <td>11.249072</td>
      <td>1008</td>
    </tr>
    <tr>
      <th>C6847</th>
      <td>24.170331</td>
      <td>4.931951</td>
      <td>4.217817</td>
      <td>5.048022</td>
      <td>0.000000</td>
      <td>11.743821</td>
      <td>7.418318</td>
      <td>4.117901</td>
      <td>1008</td>
    </tr>
    <tr>
      <th>C6845</th>
      <td>17.704774</td>
      <td>7.381424</td>
      <td>5.580256</td>
      <td>9.209765</td>
      <td>11.343971</td>
      <td>11.369700</td>
      <td>6.066916</td>
      <td>13.265252</td>
      <td>1007</td>
    </tr>
    <tr>
      <th>C6843</th>
      <td>36.751507</td>
      <td>16.630778</td>
      <td>15.052060</td>
      <td>19.661747</td>
      <td>19.661914</td>
      <td>20.375645</td>
      <td>15.319248</td>
      <td>18.352717</td>
      <td>1001</td>
    </tr>
  </tbody>
</table>
<p>168 rows × 9 columns</p>
</div>




```python
normPivAvgCov = pivAvgCov.filter(regex='A.tig')
normPivAvgCov = normPivAvgCov.apply(lambda row: scale(row, copy=True), axis=1)
normPivAvgCov.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>F_A.tig_17118</th>
      <th>F_A.tig_21321</th>
      <th>F_A.tig_21544</th>
      <th>F_A.tig_21546</th>
      <th>M_A.tig_17115</th>
      <th>M_A.tig_21322</th>
      <th>M_A.tig_21323</th>
      <th>M_A.tig_21545</th>
    </tr>
    <tr>
      <th>scaffold</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>scaffold34</th>
      <td>1.360612</td>
      <td>-0.434942</td>
      <td>0.417269</td>
      <td>0.531072</td>
      <td>-2.181267</td>
      <td>-0.421826</td>
      <td>0.752980</td>
      <td>-0.023898</td>
    </tr>
    <tr>
      <th>scaffold36</th>
      <td>1.355653</td>
      <td>-0.363497</td>
      <td>0.603710</td>
      <td>0.255252</td>
      <td>-2.156132</td>
      <td>-0.505271</td>
      <td>0.834078</td>
      <td>-0.023793</td>
    </tr>
    <tr>
      <th>C7127</th>
      <td>-0.737203</td>
      <td>0.695486</td>
      <td>-0.737203</td>
      <td>1.090714</td>
      <td>1.899816</td>
      <td>-0.737203</td>
      <td>-0.737203</td>
      <td>-0.737203</td>
    </tr>
    <tr>
      <th>C7125</th>
      <td>0.049898</td>
      <td>-0.553445</td>
      <td>0.932212</td>
      <td>-2.053928</td>
      <td>0.484745</td>
      <td>-0.227574</td>
      <td>1.514998</td>
      <td>-0.146906</td>
    </tr>
    <tr>
      <th>C7123</th>
      <td>2.339676</td>
      <td>-0.479571</td>
      <td>0.181414</td>
      <td>-0.990099</td>
      <td>0.121653</td>
      <td>0.176991</td>
      <td>-0.271759</td>
      <td>-1.078305</td>
    </tr>
  </tbody>
</table>
</div>




```python
fig, ax = plt.subplots(figsize=(10,15))
ax = sns.heatmap(normPivAvgCov, 
                 yticklabels=False,
                 ax=ax)
plt.xticks(rotation=90);
```


![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_143_0.png)



```python
g = sns.clustermap(normPivAvgCov.T,
               yticklabels=True,
               col_cluster=True,
                row_cluster=False,
#                figsize=(25, 8)
                  )

hm = g.ax_heatmap.get_position()
g.ax_heatmap.add_patch(Rectangle((19, 0.01), 5, 7.95, fill=False, edgecolor='black', lw=2))


g.ax_heatmap.set_position([hm.x0-0.1, hm.y0, hm.width*4, hm.height])
col = g.ax_col_dendrogram.get_position()
g.ax_col_dendrogram.set_position([col.x0-0.1, col.y0, col.width*4, col.height*0.5])

plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0);
```

    /home/dut/anaconda3/lib/python3.6/site-packages/matplotlib/cbook.py:136: MatplotlibDeprecationWarning: The axisbg attribute was deprecated in version 2.0. Use facecolor instead.
      warnings.warn(message, mplDeprecation, stacklevel=1)



![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_144_1.png)



```python
heatMapMaleScaffolds = ['C6889', 'C6885', 'C6883', 'C7069','C6925']
heatMapMaleScaffoldsStr = ' '.join(heatMapMaleScaffolds)
```

#### Data For Candidate Scaffolds


```python
canDf = allSampsAvgDf[allSampsAvgDf.scaffold.isin(heatMapMaleScaffolds)]
canDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>avg_cov</th>
      <th>sample_name</th>
      <th>scaffold</th>
      <th>sex</th>
      <th>length</th>
      <th>sex_name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>23506</th>
      <td>17.442269</td>
      <td>A.tig_21322</td>
      <td>C7069</td>
      <td>M</td>
      <td>1388</td>
      <td>M_A.tig_21322</td>
    </tr>
    <tr>
      <th>23507</th>
      <td>0.000000</td>
      <td>A.tig_21544</td>
      <td>C7069</td>
      <td>F</td>
      <td>1388</td>
      <td>F_A.tig_21544</td>
    </tr>
    <tr>
      <th>23508</th>
      <td>26.863238</td>
      <td>A.tig_17115</td>
      <td>C7069</td>
      <td>M</td>
      <td>1388</td>
      <td>M_A.tig_17115</td>
    </tr>
    <tr>
      <th>23509</th>
      <td>0.000000</td>
      <td>A.tig_21546</td>
      <td>C7069</td>
      <td>F</td>
      <td>1388</td>
      <td>F_A.tig_21546</td>
    </tr>
    <tr>
      <th>23510</th>
      <td>0.000000</td>
      <td>A.tig_17118</td>
      <td>C7069</td>
      <td>F</td>
      <td>1388</td>
      <td>F_A.tig_17118</td>
    </tr>
  </tbody>
</table>
</div>




```python
sns.barplot(x='scaffold', y='avg_cov', hue='sex',data = canDf)
```




    <matplotlib.axes._subplots.AxesSubplot at 0x7f30c69a45f8>




![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_148_1.png)



```python
df_list = []
for key in indvBedGraphs:
    mapped_reads = normFactorsDF[normFactorsDF['sample_name']==key[0]]['mapped_reads'].iloc[0]
    df = pd.read_csv(indvBedGraphs[key], sep='\t', names = ['scaffold', 'start', 'stop', 'coverage'])
    df = df[df.scaffold.isin(heatMapMaleScaffolds)]
    df['sample_name'] = key[0]
    df['sex'] = key[1]
    df['sex_name'] = '_'.join(key)
    df_list.append(df)
    df = None

canScaffoldsCovDf =pd.concat(df_list).reset_index(drop=True)
canScaffoldsCovDf.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>scaffold</th>
      <th>start</th>
      <th>stop</th>
      <th>coverage</th>
      <th>sample_name</th>
      <th>sex</th>
      <th>sex_name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>C6883</td>
      <td>0</td>
      <td>3</td>
      <td>2</td>
      <td>A.tig_21545</td>
      <td>M</td>
      <td>A.tig_21545_M</td>
    </tr>
    <tr>
      <th>1</th>
      <td>C6883</td>
      <td>3</td>
      <td>9</td>
      <td>3</td>
      <td>A.tig_21545</td>
      <td>M</td>
      <td>A.tig_21545_M</td>
    </tr>
    <tr>
      <th>2</th>
      <td>C6883</td>
      <td>9</td>
      <td>43</td>
      <td>4</td>
      <td>A.tig_21545</td>
      <td>M</td>
      <td>A.tig_21545_M</td>
    </tr>
    <tr>
      <th>3</th>
      <td>C6883</td>
      <td>43</td>
      <td>89</td>
      <td>5</td>
      <td>A.tig_21545</td>
      <td>M</td>
      <td>A.tig_21545_M</td>
    </tr>
    <tr>
      <th>4</th>
      <td>C6883</td>
      <td>89</td>
      <td>115</td>
      <td>6</td>
      <td>A.tig_21545</td>
      <td>M</td>
      <td>A.tig_21545_M</td>
    </tr>
  </tbody>
</table>
</div>




```python
maleScaffoldLensDF[maleScaffoldLensDF.scaffold == scaffold].length.iloc[0]
```




    1044




```python
colorsDict = dict(zip(['F','M'], ['#66c2a5', '#fc8d62'])) 
for scaffold in canScaffoldsCovDf.scaffold.unique():
    fig, axarr = plt.subplots(8, 1, figsize=(15,10))
    fig.suptitle(scaffold, fontsize=16, y=0.90)
    for i, sex_name in enumerate(sorted(canScaffoldsCovDf.sex_name.unique(), key=lambda x: x[-1])):
        data = canScaffoldsCovDf[(canScaffoldsCovDf.scaffold==scaffold) & (canScaffoldsCovDf.sex_name == sex_name) & (canScaffoldsCovDf.coverage > 0)].reset_index(drop=True)
        ax = axarr[i]
        ax.set_ylim(0,30)
        ax.set_xlim(0, maleScaffoldLensDF[maleScaffoldLensDF.scaffold == scaffold].length.iloc[0])
        ax.hlines(y=data.coverage, xmin=data.start, xmax=data.stop, color=colorsDict[sex_name[-1]], label=sex_name)
        for k,x in enumerate(zip(data.start, data.stop)):
            ax.fill_between(x=x, y1=data.coverage.iloc[k], y2=0, color=colorsDict[sex_name[-1]],alpha=1.0, lw=0)
        if i != 7:
            ax.xaxis.set_ticklabels([])
        else:
            ax.set_xlabel('Position (bp)')
        ax.legend() 
```


![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_151_0.png)



![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_151_1.png)



![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_151_2.png)



![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_151_3.png)



![png](marm_dnaseq_sex_determ_files/marm_dnaseq_sex_determ_151_4.png)



```bash
%%bash -s "$heatMapMaleScaffoldsStr"
samtools faidx ../data/male_assembly/male_assembly.scafSeq $1 > ../data/male_candidate_sequences.fasta
cat ../data/male_candidate_sequences.fasta
```

    >C6889
    TTTGCCAAACCATGTGTTTTTTTGGTCTGAACTGGACCTGGCTGGTCCAGACACAGTGGT
    TTGGGCCATGAAGAGGTAAGTGCACAAAGGAGATGCTTAAGAAATTCCAGCACTAATTCA
    GGTTCTTTAACTGCTGAGAAACTGAATTGGGATATCTTCAATTTTCATGCTTGGTCATAG
    ATTCAGCATAGCTCTTGAGTGAAAATCAAGTTTATTCTGACCATTCTTGTTCTAAGTTTT
    TTAAAAACTAAAATTTTCATACTAATTTTTTTTCCAATTTAGGGATTAGCACAAAGACAC
    ACCGTCAATGAACCCTTGCAAAAGTACGTGGAAGAAGGCTGATTGATCCATTCGTATGAT
    TTTTAGGCATTGTGAAAATCACTAAAAACCATTGAGTTCTATATCCCATGTCCCTTTTTT
    TTCATGAATAATGGGCCATGGTATAGTTGAGACAATAAGAATATATCTTGATATTTCTCA
    GACTCTCAGAAGTTTTTGATACTGTGCAGATACTTTCTGTGCAATAGTTCTGGTTCTGCC
    TTAGGGTCTGTGATCAGAAAATAATATTGATCCTTAGGGGTTCTGTAAGTTTCACTCCTA
    TCTCCTATACTATTTAAAATTTACACATATCCCTTAGGAGCAGTCATCCAGAATTTTAGA
    CTGAACAGGCAACCCTATGGAGATGACACTCAGTTCTGTCACACCTTTATGAGGTCTTTT
    ACCCTGACCTGTAGCTGGGGGAACAGCAGGCCCTCCCAATGCAAGTACTGTCGCCATTTT
    AGCTACTTAATTTATATGTTTCCATACTGAATTTATTGGGGCCTGCTCATAGACCTTTGG
    CTTGAGGCAGGATAGAAATGTAAATAAACAAGTGATAAATAAAACAAGAGCTTGGGGCTG
    AAGTGCTGAAGCAGAAATTGTAAGCTGCAATGGACTGTATGTTAACCAATATATTGCAGC
    TCAGTCTTGGTTAGATGGAATTTCTTTGCGTGCCCAGTGCTTGAGTCTAGGAGCTGGCTA
    TGACTCTGCTGTGAATGGAGTGCCTC
    >C6885
    TCTGGCCAGGAAATGTGCACTTTTACATGGATTGGCTGAAGAGGAAAAAAGAAACCCACC
    CTGCACTATCATGGAAATTTTTCTTCTCCTTGAGATGTTAAGCCTTACTCCCTTCCCTGG
    CCATTCTCCCCCTGGTTTATCTGTCAAAATCGATAAGATTCCTGACCTCCAAATTGTTTG
    TACCTCTGCATTTGCAGCTCTGAAAACAAGCAATCAAGGGAAGGGAGGGACTCCTAGGTG
    AAGAAAGCAAGCAAAGCAAGCAAGTCGGCTGAACTGAAGAAGAAGGTGAAGAAAGGAAAA
    GACATATGAAGTTCAAGCAAAGAGGGAGAGGGGAGCACGGATCAAAGAGGCCAAAAAGAG
    GAAAGAGAAAAGAGAAGCACCAATGGAGTTAGAAGGGGCTGTGCAGTAGGAAGGGAGGAA
    GTCATTCTCCAGAGTAAAGAAAAGTGAAGGTAATCCAGCAAAAAGACAGAGAAGGTCAGT
    CACTGCTCCTTGGCTGGCATCCCTTTTTACCAAGTTAGCCAATGCTTCCCTCTTGGCCCA
    ACTCAACCTTTCCCAGTGACACTTGTGCTGAACAGCACAAGGATTTCTGATCTCAGAGTG
    CTACTTAGTCTCCCCCACACCAGTGCTTAACAGTTACTGCCCTCTATAGTTGTCCCATAA
    CTGCTGAGATGGGCAAATATTGTTTCTGTACAAGAAGAAAACTAAAGAGGGGAAATCTAT
    TTCAATAAATAGATCTTGAATAACATGATAGCTTTATTAGAGGCAACAGAATTTTTGTTG
    ATTGTGTGGGTTTTGTTTTGTTTTTTATTCCTTTATTCCTCCTCCCCCCCCAAAAAAAGA
    CACCTACAAGTTCCAGCACTGGCTGAACAAGGCAACCTAAACCTTCTAATTGTGGATAGC
    ATGATAACTGCCCCAAAGGTTCTCTTTCATGTTCACATTTTTCAATAAAAATGGGCTAAA
    ATGCAGGAACGGTTTCAGTACCGTAATAATCTGATTATAGGTAAGTACTGAAGTGGAGGA
    GAAAGCTCAGAGAAAAATCAGAACT
    >C6883
    TCTGGCCAGGAAATGTGCACTTTTACATGGATTGGCTGAAGAGGAAAAAAGAAACCCACC
    CTGCACTATCATGGAAATTTTTCTTCTCCTTGAGATGTTAAGCCTTACTCCCTTCCCTGG
    CCATTCTTCCCCTGGTTTATCTGTCAAAATCGATAAGATTCCTGACCTCCAAATTGTTTG
    TACCTCTGCATTTGCAGCTCTGAAAACAAGCAATCAAGGGAAGGGAGGGACTCCTAGGTG
    AAGAAAGCAAGCAAAGCAAGCAAGTCGGCTGAACTGAAGAAGAAGGTGAAGAAAGGAAAA
    GACATATGAAGTTCAAGCAAAGAGGGAGGGGGAGCACGGATCAAAGAGGCCAAAAAGAGG
    AAAGAGAAAAGAGAAGCACCAATGGAGTTAGAAGGGGCTGTGCAGTAGGAAGGGAGGAAG
    TCATTCTTCAGAGTAAAGAAAAGTGAAGGTAATCCAGCAAAAAGACAGAGAAGGTCAGTT
    ACTGCTCCTTGGCTGGCATCCCTTTTTACCAAGTTAGCCAATGCTTCCCTCTTGGCCCAA
    CTCAACCTTTCCCAGTGACACTTGTGCTGAACAGCACAAGGATTTCTGATCTCAGAGTGC
    TACTTAGTCTCCCCCACACCAGTGCTTAACAGTAACTGCCCTCTATAGTTGTCCCATAAC
    TGCTGAGATGGGCAAATATTGTTTCTGTACAAGAAGAAAACTAAAGAGGGGAAATCTATT
    TCAATAAATAGATCTTGAATAACATGATAGCTTTATTAGAGGCAACAGAATTTTTGTTGA
    TTGTGTGGGTTTTGTTTTGTTTTTTATTCCTTTATTCCTCCTCCCCCCCAAAAAAAAGAC
    ACCTACAAGTTCCAGCACTGGCTGAACAAGGCAACCTAAACCTTCTAATTGTGGATAGCA
    TGATAACTGCCCCAAAAGTTCTCTTTCATGTTCACATTTTTCAATAAAAATGGGCTAAAA
    TGCAGGAACGGTTTCAGTACCGTAATAATCTGATTATAGGTAAGTACTGAAGTGGAGGAG
    AAAGCTCAGAGAAAAATCAGAACT
    >C7069
    AGTCCACCAGCCGCCACTGACTACAAGCCACCTGAAAAGAAGGAGGAAAGGGTTGAAAAA
    AAAAGGGCATACTGCCAACACCTTTCAAACTGAAAACAGATTATAATCTGATCCAGATAG
    AGGAAAATCCATTTGGAGAGAGTACAGTCGAGGGGAAAAAAATCCAGCTTAGTCTTTAAA
    AGGTGTAACTGTGGCATTTTAAGTTTTCTTGATATCAGCAATATTAGGAATGGCTGCACC
    TAGCTGATCAGATCTGGATGGGTTGATTTGATTATTGCTTAAATTGTCATCAAACCCTGG
    ATTCACCATAGTCAAGTAAAAGCTGTAAAAACACTGGATACCATAGAGATGTTTGTATAT
    TTTTATTTAAAATTTAATGCCACTTTGTACAGCAAACTGTTCACTTTGATATATGATTGT
    GCCTTACAAACGAAAATAAGTTTTTTCACTTTCCTTAACATTTTAAATTTGAATTTCTTA
    TGCCTAGAATGATATGCAAATGACACAACTACAATGAAATCAAAAGTGTTACAAAATGTC
    TACAAAGGAATAATCCCATGAATGTCTCAGATAACAGTTGTATTAATTATGATGGAGCAA
    CAGGTAAAGAGACAGTGAAGCGAAAGGACATTTTCACCTCAAAGGTTGGTGAAGAAACTG
    CATCATAGCACATCTCCAATTTCAAGGGTTAGTAATGAAGTGATGCCAGTTTTTGCAGCT
    GTAGGAAAGTATGCTCATGCAGATCTAGGCAACTATACTGGTGATCCCAGGTTTTCTAAT
    TGAAGAAGAAAATCACAGCTGGGGAAAATAATGGCAAAGTATGGCACAAGAGGAGAAGTT
    GTTCCCTTAAAGAAATTACTTTGCCATTGCCTCCCCCTGCAAGTTCACAAAATTAATAGT
    ATACTCTTGAGTTTATTTTGAAAACAAAGAAAGACCTAGTGGCAGGTTTATAACTGTGTC
    CGAAACGTCAGGTTGCTTCCTGTCAATGACCGAGCACAGACATCATCAATCCATTCATTT
    TGTGCAACTGTGAAATGATTTTTCCAATCTCCAACAATGCCTGCAAATGAGAGAGACCAT
    CAAGTCTTGTTTTTAATATGATATTACGCAAGCAGTTACATTAATTGGTAATAGGTTTTT
    TAATTGATTAATTCCTCAGCAATTCTACAATTCAGCATTAATACATATGTTATCAGCACA
    CGTGATTTTTTTTTTAAAAAATCTGATTTAATTGAATTTGTTTGACAGCAGTTTATACAT
    GAAGAAGAAGAACAGAAGGACTAGGCTTTTATATCCACTTTTCTCATCCTTAAGGAATCT
    CAGGACACCCTGTGTGGGGGGGGCTGCAAAAGCTCAGAAGGAGAGCTATGATTACCCAAG
    GCCACCCA
    >C6925
    ATAAATAATGGGAACATTGGATACTGAGTTCAAATTCACCTCGAATTTTAGTGCCTATGG
    AAATGCACCCTAGGTTTCTTTCTTATGATTTCTTCATCAACAGGGACTACCTGGATGCAG
    GAAATTGTAGACATGATCCAACATGGGGGAGACCCCCAGCAATGTGCTCGGGCTCCTATC
    TATGAAAGGAATCCATTCATAGAGCTGTGTCCTCCAAAACCTATTCGTACAGGTGAGAAA
    CAAGTAGTACAGAGTTCCAGTGAAGTCATGCATGAAAACAGAGTTACTTAAGATGAATCT
    GCATATTCCTCTCTGTTTCAGGTTTTGAAGAAGCAGAGGTAATGCCTTCCCCACGCACAC
    TCAAATCACACCTCCCTGTCCATCTTTTACCACCTTCCTTCTGGGAACAGAATTGTAAGG
    TCAGAAGAATAAGTGGGGATCAAGTTGAAGACAGGACTGATGTCCAGACTGATTCTCACT
    TCATTGTGAACCATCAACTTATTAATAATGTGAAAAACAATGGTCTTTCAAAATTAATTG
    TGCCACAGACTTGGCAATGAAATTTTCCCTCCGCAAAACACAGTCAGAAATTGATACCTT
    ATAAAAATGGTGCATTACTGAAAATAATGAAAAGAGGATAATACCAGGGAAAAAATTATT
    ACAGAAGTGAAAACTTTTTAAAAAGAGGTACAAAGTGTTTAGATATATGACCTGATAAAC
    CAACATTTCACCTAGTCTAAACATCTTCAAACTAAAGACTTAATTAGAATTACTGTTACT
    GCTTGTTCTTTTGCAGAACTTTGTAATACAATCATATCATAATCCTAGTAACGCTACTGA
    TAACAGTACAATTATGGTTTTGGATCAGGGGTTGTGTTGTGTACATACGTGCTTGCAAGG
    TCATGGCTGGGGAGGTGGAAGACCATTCAATTTGTGGGCAAAGGGGGGACCCTACAGATG
    TAGGTCCTTTTTGATTCACCTGCAAAAGTTCAAATTAATAGGAATTCATACAAATTGCAA
    ATATACTGATTATATTTTGCCGCATTTCTGTAGATCATTTATGTAGCCAGGAATGTCAAG
    GACACTGCA


#### Designing Primers


```python
%load_ext autoreload

%autoreload 2
```

    The autoreload extension is already loaded. To reload it, use:
      %reload_ext autoreload



```python
import fasta_classes as fa
```


```python
targetSeqs = fa.Fasta_file('../data/male_candidate_sequences.fasta')
targeSeqTups = list(targetSeqs.fasta_tuples)
```


```python
targeSeqTups[0]
```




    ('>C6889',
     'TTTGCCAAACCATGTGTTTTTTTGGTCTGAACTGGACCTGGCTGGTCCAGACACAGTGGTTTGGGCCATGAAGAGGTAAGTGCACAAAGGAGATGCTTAAGAAATTCCAGCACTAATTCAGGTTCTTTAACTGCTGAGAAACTGAATTGGGATATCTTCAATTTTCATGCTTGGTCATAGATTCAGCATAGCTCTTGAGTGAAAATCAAGTTTATTCTGACCATTCTTGTTCTAAGTTTTTTAAAAACTAAAATTTTCATACTAATTTTTTTTCCAATTTAGGGATTAGCACAAAGACACACCGTCAATGAACCCTTGCAAAAGTACGTGGAAGAAGGCTGATTGATCCATTCGTATGATTTTTAGGCATTGTGAAAATCACTAAAAACCATTGAGTTCTATATCCCATGTCCCTTTTTTTTCATGAATAATGGGCCATGGTATAGTTGAGACAATAAGAATATATCTTGATATTTCTCAGACTCTCAGAAGTTTTTGATACTGTGCAGATACTTTCTGTGCAATAGTTCTGGTTCTGCCTTAGGGTCTGTGATCAGAAAATAATATTGATCCTTAGGGGTTCTGTAAGTTTCACTCCTATCTCCTATACTATTTAAAATTTACACATATCCCTTAGGAGCAGTCATCCAGAATTTTAGACTGAACAGGCAACCCTATGGAGATGACACTCAGTTCTGTCACACCTTTATGAGGTCTTTTACCCTGACCTGTAGCTGGGGGAACAGCAGGCCCTCCCAATGCAAGTACTGTCGCCATTTTAGCTACTTAATTTATATGTTTCCATACTGAATTTATTGGGGCCTGCTCATAGACCTTTGGCTTGAGGCAGGATAGAAATGTAAATAAACAAGTGATAAATAAAACAAGAGCTTGGGGCTGAAGTGCTGAAGCAGAAATTGTAAGCTGCAATGGACTGTATGTTAACCAATATATTGCAGCTCAGTCTTGGTTAGATGGAATTTCTTTGCGTGCCCAGTGCTTGAGTCTAGGAGCTGGCTATGACTCTGCTGTGAATGGAGTGCCTC')




```python
results_dfs = []
for seq in targeSeqTups:
    df = fa.design_primers(seq, fragment_size=[300,500],exlude_seqs='../../dovetail_genome_delivery/data/fastas/lizard_23Jun2015_piz6a.upper.fasta')
    results_dfs.append(df)

allPrimerDf= pd.concat(results_dfs).reset_index(drop=True)
allPrimerDf
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>SEQ_NAME</th>
      <th>PRIMER_SET</th>
      <th>PRIMER_LEFT_SEQUENCE</th>
      <th>PRIMER_RIGHT_SEQUENCE</th>
      <th>PRIMER_PAIR_PRODUCT_SIZE</th>
      <th>PRIMER_LEFT_COORD</th>
      <th>PRIMER_RIGHT_COORD</th>
      <th>PRIMER_LEFT_GC_PERCENT</th>
      <th>PRIMER_RIGHT_GC_PERCENT</th>
      <th>PRIMER_LEFT_TM</th>
      <th>PRIMER_RIGHT_TM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>C6889</td>
      <td>0</td>
      <td>ACTGAACAGGCAACCCTATGG</td>
      <td>ACTGGGCACGCAAAGAAATTC</td>
      <td>339</td>
      <td>(659, 21)</td>
      <td>(997, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.995255</td>
      <td>60.000274</td>
    </tr>
    <tr>
      <th>1</th>
      <td>C6889</td>
      <td>1</td>
      <td>TCCCTTAGGAGCAGTCATCCA</td>
      <td>ACTGGGCACGCAAAGAAATTC</td>
      <td>369</td>
      <td>(629, 21)</td>
      <td>(997, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.991673</td>
      <td>60.000274</td>
    </tr>
    <tr>
      <th>2</th>
      <td>C6889</td>
      <td>2</td>
      <td>CTGGTTCTGCCTTAGGGTCTG</td>
      <td>ACTGGGCACGCAAAGAAATTC</td>
      <td>469</td>
      <td>(529, 21)</td>
      <td>(997, 21)</td>
      <td>57.142857</td>
      <td>47.619048</td>
      <td>60.065854</td>
      <td>60.000274</td>
    </tr>
    <tr>
      <th>3</th>
      <td>C6885</td>
      <td>0</td>
      <td>AAGCAATCAAGGGAAGGGAGG</td>
      <td>AAAGGTTGAGTTGGGCCAAGA</td>
      <td>344</td>
      <td>(207, 21)</td>
      <td>(550, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.993600</td>
      <td>60.063266</td>
    </tr>
    <tr>
      <th>4</th>
      <td>C6885</td>
      <td>1</td>
      <td>AAGCAATCAAGGGAAGGGAGG</td>
      <td>AGCACAAGTGTCACTGGGAAA</td>
      <td>362</td>
      <td>(207, 21)</td>
      <td>(568, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.993600</td>
      <td>60.064758</td>
    </tr>
    <tr>
      <th>5</th>
      <td>C6885</td>
      <td>2</td>
      <td>GAAACCCACCCTGCACTATCA</td>
      <td>GGAGCAGTGACTGACCTTCTC</td>
      <td>439</td>
      <td>(50, 21)</td>
      <td>(488, 21)</td>
      <td>52.380952</td>
      <td>57.142857</td>
      <td>59.995723</td>
      <td>60.067213</td>
    </tr>
    <tr>
      <th>6</th>
      <td>C6883</td>
      <td>0</td>
      <td>GAAACCCACCCTGCACTATCA</td>
      <td>CAGCCAAGGAGCAGTAACTGA</td>
      <td>445</td>
      <td>(50, 21)</td>
      <td>(494, 21)</td>
      <td>52.380952</td>
      <td>52.380952</td>
      <td>59.995723</td>
      <td>59.997570</td>
    </tr>
    <tr>
      <th>7</th>
      <td>C6883</td>
      <td>1</td>
      <td>GAAACCCACCCTGCACTATCA</td>
      <td>AGCCAAGGAGCAGTAACTGAC</td>
      <td>444</td>
      <td>(50, 21)</td>
      <td>(493, 21)</td>
      <td>52.380952</td>
      <td>52.380952</td>
      <td>59.995723</td>
      <td>59.997456</td>
    </tr>
    <tr>
      <th>8</th>
      <td>C6883</td>
      <td>2</td>
      <td>GGCCATTCTTCCCCTGGTTTA</td>
      <td>CAGCCAAGGAGCAGTAACTGA</td>
      <td>377</td>
      <td>(118, 21)</td>
      <td>(494, 21)</td>
      <td>52.380952</td>
      <td>52.380952</td>
      <td>59.993361</td>
      <td>59.997570</td>
    </tr>
    <tr>
      <th>9</th>
      <td>C7069</td>
      <td>0</td>
      <td>GCTGCACCTAGCTGATCAGAT</td>
      <td>CAGCTGCAAAAACTGGCATCA</td>
      <td>489</td>
      <td>(232, 21)</td>
      <td>(720, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.929555</td>
      <td>60.269883</td>
    </tr>
    <tr>
      <th>10</th>
      <td>C7069</td>
      <td>1</td>
      <td>GCTGCACCTAGCTGATCAGAT</td>
      <td>GCTGCAAAAACTGGCATCACT</td>
      <td>487</td>
      <td>(232, 21)</td>
      <td>(718, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.929555</td>
      <td>60.269936</td>
    </tr>
    <tr>
      <th>11</th>
      <td>C7069</td>
      <td>2</td>
      <td>GCTGCACCTAGCTGATCAGAT</td>
      <td>AGCTGCAAAAACTGGCATCAC</td>
      <td>488</td>
      <td>(232, 21)</td>
      <td>(719, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.929555</td>
      <td>60.269936</td>
    </tr>
    <tr>
      <th>12</th>
      <td>C6925</td>
      <td>0</td>
      <td>TGCTCGGGCTCCTATCTATGA</td>
      <td>CTGTGTTTTGCGGAGGGAAAA</td>
      <td>419</td>
      <td>(164, 21)</td>
      <td>(582, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>59.924464</td>
      <td>59.592742</td>
    </tr>
    <tr>
      <th>13</th>
      <td>C6925</td>
      <td>1</td>
      <td>GAGCTGTGTCCTCCAAAACCT</td>
      <td>CTGTGTTTTGCGGAGGGAAAA</td>
      <td>382</td>
      <td>(201, 21)</td>
      <td>(582, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>60.202818</td>
      <td>59.592742</td>
    </tr>
    <tr>
      <th>14</th>
      <td>C6925</td>
      <td>2</td>
      <td>AGAGCTGTGTCCTCCAAAACC</td>
      <td>CTGTGTTTTGCGGAGGGAAAA</td>
      <td>383</td>
      <td>(200, 21)</td>
      <td>(582, 21)</td>
      <td>52.380952</td>
      <td>47.619048</td>
      <td>60.202818</td>
      <td>59.592742</td>
    </tr>
  </tbody>
</table>
</div>




```python
temp = allPrimerDf[["SEQ_NAME", "PRIMER_LEFT_SEQUENCE", "PRIMER_LEFT_GC_PERCENT", 'PRIMER_LEFT_TM']].copy()

temp['sequence'] = temp["PRIMER_LEFT_SEQUENCE"]
temp['description/comments'] = temp.apply(lambda x: 'Left (top) Probe for {}'.format(x['SEQ_NAME']), axis=1)
temp['length (nt)'] = temp['sequence'].apply(lambda x: len(x))
temp['requested by'] = 'dut'
temp['%GC'] = temp['PRIMER_LEFT_GC_PERCENT']
temp['strand'] = 'T'
temp['Tm'] = temp['PRIMER_LEFT_TM']
temp = temp[['sequence', 'description/comments', 'length (nt)', 'requested by', '%GC','strand', 'Tm']]


temp2 = allPrimerDf[["SEQ_NAME", "PRIMER_RIGHT_SEQUENCE", "PRIMER_RIGHT_GC_PERCENT", "PRIMER_RIGHT_TM"]].copy()


temp2['sequence'] = temp2["PRIMER_RIGHT_SEQUENCE"]
temp2['description/comments'] = temp2.apply(lambda x: 'Right (bottom) Probe for {}'.format(x['SEQ_NAME']), axis=1)
temp2['length (nt)'] = temp2['sequence'].apply(lambda x: len(x))
temp2['requested by'] = 'dut'
temp2['%GC'] = temp2['PRIMER_RIGHT_GC_PERCENT']
temp2['strand'] = 'B'
temp2['Tm'] = temp2['PRIMER_RIGHT_TM']
temp2 = temp2[['sequence', 'description/comments', 'length (nt)', 'requested by', '%GC','strand', 'Tm']]

order_df = pd.concat((temp, temp2)).drop_duplicates().reset_index(drop=True)
order_df = order_df.reset_index()

order_df['name'] = order_df.apply(lambda x: 'Holi{}{}'.format(int(x['index']) + 1252, x['strand']), axis=1)
order_df = order_df.drop('index',axis=1)
order_df['date'] = datetime.datetime.today().strftime('%-m/%-d/%Y')
order_df['Distance to start of primer (nt)'] = np.nan
order_df['Page Reference'] = os.getcwd() + '/marm_dnaseq_sex_determ.ipynb'
order_df['Additional Comments/Results'] = 'Primers for Male Specific Sequences in Aspidoscelis marmoratus'
order_df['%GC'] = round(order_df['%GC'],3)
order_df['Tm'] = round(order_df['Tm'],3)
order_df = order_df[['name', 'sequence', 'length (nt)','date', 'Distance to start of primer (nt)','description/comments', 'Page Reference', 'Tm', 'Additional Comments/Results', '%GC']]

order_df
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>sequence</th>
      <th>length (nt)</th>
      <th>date</th>
      <th>Distance to start of primer (nt)</th>
      <th>description/comments</th>
      <th>Page Reference</th>
      <th>Tm</th>
      <th>Additional Comments/Results</th>
      <th>%GC</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Holi1252T</td>
      <td>ACTGAACAGGCAACCCTATGG</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6889</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.995</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Holi1253T</td>
      <td>TCCCTTAGGAGCAGTCATCCA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6889</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.992</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Holi1254T</td>
      <td>CTGGTTCTGCCTTAGGGTCTG</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6889</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.066</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>57.143</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Holi1255T</td>
      <td>AAGCAATCAAGGGAAGGGAGG</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6885</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.994</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Holi1256T</td>
      <td>GAAACCCACCCTGCACTATCA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6885</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.996</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Holi1257T</td>
      <td>GAAACCCACCCTGCACTATCA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6883</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.996</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Holi1258T</td>
      <td>GGCCATTCTTCCCCTGGTTTA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6883</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.993</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Holi1259T</td>
      <td>GCTGCACCTAGCTGATCAGAT</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C7069</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.930</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Holi1260T</td>
      <td>TGCTCGGGCTCCTATCTATGA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6925</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.924</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Holi1261T</td>
      <td>GAGCTGTGTCCTCCAAAACCT</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6925</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.203</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>10</th>
      <td>Holi1262T</td>
      <td>AGAGCTGTGTCCTCCAAAACC</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Left (top) Probe for C6925</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.203</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>11</th>
      <td>Holi1263B</td>
      <td>ACTGGGCACGCAAAGAAATTC</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C6889</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.000</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>47.619</td>
    </tr>
    <tr>
      <th>12</th>
      <td>Holi1264B</td>
      <td>AAAGGTTGAGTTGGGCCAAGA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C6885</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.063</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>47.619</td>
    </tr>
    <tr>
      <th>13</th>
      <td>Holi1265B</td>
      <td>AGCACAAGTGTCACTGGGAAA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C6885</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.065</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>47.619</td>
    </tr>
    <tr>
      <th>14</th>
      <td>Holi1266B</td>
      <td>GGAGCAGTGACTGACCTTCTC</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C6885</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.067</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>57.143</td>
    </tr>
    <tr>
      <th>15</th>
      <td>Holi1267B</td>
      <td>CAGCCAAGGAGCAGTAACTGA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C6883</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.998</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>16</th>
      <td>Holi1268B</td>
      <td>AGCCAAGGAGCAGTAACTGAC</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C6883</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.997</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>52.381</td>
    </tr>
    <tr>
      <th>17</th>
      <td>Holi1269B</td>
      <td>CAGCTGCAAAAACTGGCATCA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C7069</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.270</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>47.619</td>
    </tr>
    <tr>
      <th>18</th>
      <td>Holi1270B</td>
      <td>GCTGCAAAAACTGGCATCACT</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C7069</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.270</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>47.619</td>
    </tr>
    <tr>
      <th>19</th>
      <td>Holi1271B</td>
      <td>AGCTGCAAAAACTGGCATCAC</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C7069</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>60.270</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>47.619</td>
    </tr>
    <tr>
      <th>20</th>
      <td>Holi1272B</td>
      <td>CTGTGTTTTGCGGAGGGAAAA</td>
      <td>21</td>
      <td>10/25/2017</td>
      <td>NaN</td>
      <td>Right (bottom) Probe for C6925</td>
      <td>/n/projects/dut/a_marmorata/dnaseq_sex_determ/...</td>
      <td>59.593</td>
      <td>Primers for Male Specific Sequences in Aspidos...</td>
      <td>47.619</td>
    </tr>
  </tbody>
</table>
</div>




```python
order_df.to_excel('../data/male_marker_primer_order_df.xls')
allPrimerDf.to_excel('../data/male_marker_primer_sets.xls')
```

## Assembling Unmapped Male Reads 67-mer


```python
male63Config ="""
maximal read length
max_rd_len=250
[LIB]
#average insert size
avg_ins=538
#if sequence needs to be reversed, 1 for forward reverse
reverse_seq=1
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=250
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#bam file for single or paired reads, reads 1 in paired reads file should always be followed by reads 2
#NOTE: If a read in bam file fails platform/vendor quality checks(the flag field 0x0200 is set), itself and it's paired read would be ignored.
b=/n/projects/dut/a_marmorata/dnaseq_sex_determ/data/malelib/male_unmapped.bam
"""
with open('../bin/male_63_Assembly.config','w') as fw:
    fw.write(male63Config)
print(male63Config)
```

    
    maximal read length
    max_rd_len=250
    [LIB]
    #average insert size
    avg_ins=538
    #if sequence needs to be reversed, 1 for forward reverse
    reverse_seq=1
    #in which part(s) the reads are used
    asm_flags=3
    #use only first 100 bps of each read
    rd_len_cutoff=250
    #in which order the reads are used while scaffolding
    rank=1
    # cutoff of pair number for a reliable connection (at least 3 for short insert size)
    pair_num_cutoff=3
    #minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
    map_len=32
    #bam file for single or paired reads, reads 1 in paired reads file should always be followed by reads 2
    #NOTE: If a read in bam file fails platform/vendor quality checks(the flag field 0x0200 is set), itself and it's paired read would be ignored.
    b=/n/projects/dut/a_marmorata/dnaseq_sex_determ/data/malelib/male_unmapped.bam
    



```bash
%%bash
cd ../bin/SOAPdenovo2-bin-LINUX-generic-r240
#nohup ./SOAPdenovo-63mer all -s ../male_63_Assembly.config -K 63 -R -p 20 -o male_63_assembly 1>ass.log 2>ass.err &
## after the assembly finished...
#mkdir ../../data/male_63_assembly
#mv male_63_assembly* ../../data/male_63_assembly/
#mv ass* ../../data/male_63_assembly/
```
