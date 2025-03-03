# CS4150 Data Science

## Activity 1: Data Exploration
### Task:

Write a python script that accepts the given data set and finds the following:

1. Number of genomic windows
2. Number of NPs
3. On average, how many windows are present in an NP?
4. What is the smallest number of windows present in any NP? The largest?
5. On average, what is the number of NPs in which a window is detected? The smallest? The largest?

## Activity 2: Radial Position and Compaction
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

## Activity 3: Nuclear Organization

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

## Activity 4: NP Similarity

### Task:

Knowing the coordinate location of the Hist1 region on mouse chromosome 13, we will study this section. To accomplish this,
we must use our extracted Hist1 region from the previous activity as well as all of their relevant NPs. Once this step has
been completed, we will then compute the Jaccard index (aka the Jaccard similarity coefficient) for relevant NPs.
Next, the computed values will be stored in an $NP$ x $NP$ matrix; containing the similarities for all pairs of relevant NPs.
[Relevant link](https://en.wikipedia.org/wiki/Jaccard_index#Similarity_of_asymmetric_binary_attributes)

## Activity 5: NP distance & heatmaps

### Task:

Continuing our study of the Hist1 region, in this activity, I used the previous activity's computed values (containing similarity scores)
to compute Jaccard distance between the relevant NPs. This result is stored in a new matrix. Finally, both similarity and distance will
be represented as a heatmap. [Relevant link](https://en.wikipedia.org/wiki/Jaccard_index#Similarity_of_asymmetric_binary_attributes)

## Activity 6: K-Means Clustering

### Task:

In this activity, we normalize and clean the data, removing any NP's which contain no window detections. The denominator in the Jaccard equation
is updated to: $$(min{|A|, |B|})$$ Where $A$ represents the total number of attributes where the attribute of A is 1. 
<br>
Where $B$ represents the total number of attributes where the attribute of B is 1. 
<br></br>
Let the numerator in the Jaccard equation = $M_{11}$, the total number of attributes where the
attribute of $A$ is 1 and the attribute of $B$ is 1.

Use k-means clustering to _form 3 clusters_ from the NPs that captured the Hist1 region.

<ol>
    <li>
        Randomly Select k=3 distinct data points. These are the initial clusters.
    </li>
    <br>
    <li>
        Measure the distance between each point and each of the 'k' clusters using the normalization procedure described above.
    </li>
    <br>
    <li>
        Assign each point to the nearest cluster.
    </li>
</ol>

## Activity 7: K-Means Improved

### Task:

<ol>
    <li>
        Using the previous activity, find the center of each cluster. For each cluster, find the NP whose average dissimilarity to all the objects
        in the cluster is minimal. These are the centers of the new clusters.
    </li>
    <li>
        Repeat until the clusters no longer change.
    </li>
    <li>
        Assess the quality of the clustering by adding up the variation within each cluster.
    </li>
</ol>

> ### **Note:** <br>
> The NP selected by this method is called the medoid (see
https://en.wikipedia.org/wiki/Medoid).
“In contrast to the k-means algorithm, k-medoids clustering chooses actual data points
as centers (medoids or exemplars), and thereby allows for greater interpretability of the
cluster centers than in k-means, where the center of a cluster is not necessarily one of
the data points (it is the average between the points in the cluster).”

## Activity 8: Feature Selection

### Task:

Extend K-Means to accomplish:

<ol>
    <li>
        Perform K-Means clustering repeatedly, using different random starting points in each iteration. Keep track of the clusters and their quality levels (compute intra-cluster similarity and intra-cluster distance).
    </li>
    <li>
        Pick the best set of clusters (based on intra-cluster metrics from part A).
    </li>
    <li>
        For the best set of clusters, visualize EACH cluster as a heatmap, wherein:
        <ol>
            <li>
                Rows = NPs
            </li>
            <li>
                Columns = genomic windows
            </li>
            <li>
                Cells = values from the segregation table
            </li>
        </ol>
    </li>
</ol>

## Activity 9: Feature Selection Part 2

### Task:

Characterize the best k-clusters in terms of genomic and epigenomic features:

#### Features:
<ol>
    <li>
        Histone Genes
    </li>
    <li>
        Lamina-associated domains (LADs)
    </li>
</ol>

Extend the program to process the feature table and determine if there is a pattern of features in each cluster.

For each cluster, perform the following:

<ol>
    <li>
        Calculate the percentage of windows in each NP that contain histone genes
    </li>
    <li>
        Show the result of step 1 in a set of three boxplots wherein:
        <ol>
            <li>
                The x-axis represents the clusters
            </li>
            <li>
                The y-axis represents the percentage of windows in an NP that contain histone genes
            </li>
        </ol>
    </li>
    <li>
        Calculate the percentage of windows in each NP that contain LADs
    </li>
    <li>
        Show the result of step 3 in a set of three boxplots following the same rules as step 2.
    </li>
</ol>

## Activity 10: Feature Selection Part 3

### Task:

Extend the program to characterize each cluster by radial position. Consider the radial positions 
of the NPs in the best set of clusters. Is there a pattern?

For each cluster:

<ol>
    <li>
        Consider the previously calculated radial position for each NP, which is an integer between 
        1 and 5. 
        <br> 1 - Strongly apical 
        <br> 2 - Somewhat apical 
        <br> 3 - Neither apical nor equatorial
        <br> 4 - Somewhat equatorial
        <br> 5 - Strongly equatorial
    </li>
    <li>
        Calculate the percentage of NPs in each radial position category.
    </li>
    <li>
        Show the result of step 2 in a bar graph:
        <ol>
            <li>
                Show five bars - one bar for each radial position category.
            </li>
            <li>
                The height of a bar is the number of NPs classified as being in a particular radial position category.
            </li>
        </ol>
    </li>
</ol>