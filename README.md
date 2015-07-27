# hi-c_helper_scripts

A set of python scripts to deal with hi-c interactions matrices, and make inputs for circos.

##Usage:

In the command line:

1. **cross_matrix.py**: uses 2 hi-c matrixes and makes a table with the interactions in comum and not in common - requires itertools
<pre><code>python cross_matrix.py [cromossome] [start] [end] [matrix1] [matrix2] [outfile]
</code></pre>


2. **getgenes.py**:This script uses the biomart server to search for genes in the interval of the genome defined by the user. - requires biomart for python
<pre><code>python getgenes.py [chromossome] [start_interval] [end_interval] [out_exon/intro] [out_label1] [out_label2]
</code></pre>

3. **merge_hic_replicates.py**: merge the hi-c replicates. only the entries present in the two replicates are present in the final set.
<pre><code>python merge_hic_replicates.py [chromossome] [replicate1.tab] [replicate2.tab] [output]
</code></pre>

4. **sam2interactions.py**: This script uses the sam files from alignments of 3 diferent runs to construct a interaction matrix. Requires numpy and python collections
Part 1:
<pre><code>python sam2interactions.py 1 [chromossome of interst] [start of the interval] [end of the interval] [samfile1] [samfile2] [samfile3] [filtred1.sam] [filtered2.sam] [filtered3.sam]
</code></pre>

Part 2:
<pre><code>python sam2interactions.py 2 [merged.sam] [chromossome of interest] [size of the windows or interval] [intracromossome_interactions.matrix] [intercromossome_interactions.matrix]
</code></pre>

5. **script_circos_color_by_number_interactions.py**: Uses hi-c matrices and transform them into the input required by circos.
<pre><code>python script_circos.py [infile] [chromossome] [start] [end] [interval] [outfile]
</code></pre>

5. **subinteractions.py**: searches for the equal target zones of two diferent intervals, and makes a matrix with only that interactions.

If you want all the targets in common, use:
<pre><code>python subinteractions.py all [matrix] [chr1] [start1] [end1] [chr2] [start2] [end2] [outfile]
</code></pre>

#if you want the targets in genes, use:
<pre><code>python subinteractions.py genes [matrix] [chr1] [start1] [end1] [chr2] [start2] [end2] [genes_label_file] [outfile]
</code></pre>


##License:

GPLv2
