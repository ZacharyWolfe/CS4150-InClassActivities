# CS4150 Data Science

## Activity 1:
### Task:

Write a python script that accepts the given data set and finds the following:

1. Number of genomic windows
2. Number of NPs
3. On average, how many windows are present in an NP?
4. What is the smallest number of windows present in any NP? The largest?
5. On average, what is the number of NPs in which a window is detected? The smallest? The largest?

## Activity 2:
### Task:

### Evaluate data quality:
Identify NPs (i.e., samples) that have an unusually high detection frequency (the detection
frequency of an NP is the percentage of genomic windows detected by the NP).

### Estimate radial position of each NP:

The [GAM](https://pmc.ncbi.nlm.nih.gov/articles/PMC5366070/#!po=5.27638) paper describes how to estimate the radial position of an NP. Use these insights to estimate the radial position for each NP on a scale of 1-5:
1) Strongly apical
2) Somewhat apical
3) Neither apical nor equatorial
4) Somewhat equatorial
5) Strongly equatorial

### Estimate compaction of each genomic window:

The [GAM](https://pmc.ncbi.nlm.nih.gov/articles/PMC5366070/#!po=5.27638) paper describes how to estimate the degree of compaction of a genomic window. Use these insights to assign each window a compaction rating between 1-10 (10 is the most condensed; 1 is the least condensed).

## Activity 3:

### Task:

We studied the Histone Locus Body (HLB), which resides in the nucleus and plays a vital role
in the production of histone genes. The HLB is formed during the cell cycle when histone genes
are needed. At other times, the HLB is deconstructed. The HLB is formed by the coalescing of
portions of the genome that contain the histone genes (which are contained on chromosomes
13, 3, 11, and 6 of the mouse). We studied study the largest portion of the HLB – the Hist1 region.

> Note: Hist1 is located on mouse chromosome 13 at the following coordinates:
>> __Start:__ 21.7 Mb
> <br>__Stop:__ 24.1 Mb

1) <br>a) Extract the Hist1 region from the segregation table.
   <br>b) Extract relevant NPs for Hist1 (NPs that contain at least one window in the region of interest).
<br><br>
2) Calculate basic statistics for the Hist1 region:
    <br>a) Number of genomic windows
    <br>b) Number of NPs
    <br>c) On average, how many windows are present in an NP?
    <br>d) What is the smallest number of windows present in any NP? The largest?
    <br>e) On average, what is the number of NPs in which a window is detected? The smallest? The largest?
<br><br>
3) What are the most common radial positions of the Hist1 region? (Based on the NPs that captured the region.)
<br><br>
4) What are the typical compactions of the windows within the Hist1 region?

## Activity 4:

### Task:

