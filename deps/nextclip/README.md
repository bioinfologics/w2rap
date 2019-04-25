NextClip
========

Nextera Long Mate Pair analysis and processing tool. See manual for installation and instructions.

The 'examples' directory contains example reports and a partial set of reads for a (relatively low quality) Nextera library. See the ReadMe file inside the examples directory for details of running the examples.

To hear about updates to NextClip, please follow @richardmleggett on Twitter: https://twitter.com/richardmleggett

## Changes in v1.3

**Bug fixes**
- Improved ability to spot junction adaptor due to bug in previous logic.
- Removed some crashes that can occur if some of the stats counts are zero (ie. dividing zero by number of reads creates error).
- Corrected LSF submission to request correct amount of memory for NextClip (and provided option which will report memory requirements before running NextClip).
- The min_length parameter was 1 base out (i.e. specifying 25 gave reads of 24) - this is fixed, meaning that yields will appear to go down slightly. You can always decrease min_length to 24 if you want to recover this.
- Fixed PCR duplicate stats being wrong for counts > 99.
- Fixed help wording.

**Improvements**
- NextClip outputs number of bases written for each category and the pipeline PDF report shows this as well as approximate genome coverage.
- Ability to pass parameters through to the pipeline (e.g. minimum width, trim ends etc.).
- You can specify the pipeline stage to start at, so no need to redo all analysis in the event of a fail.
- The scripts directory is found automatically - no need to set it during installation.
- Moved PCR duplicate kmer offset to avoid lower quality bases at start of read. Should provide slightly more accurate estimates of PCR duplication.
- Queue name can be specified when running with LSF.
- Reduced size of hash table used to store PCR duplicates.
