# VCF-comparator
Provide a tool to compare positions presents inside variant call format (VCF files).
This tool provide two way of analyzing VCF files :
- By comparing files in there entirety (Files comparison).
- By summarizing which variants of each file are the most shared (Variants summarization) (default mod).

You can use this tools to solve a number of issues related to the replicates of an experience :
- Identification of aberrant replicates. (Files comparison)
- Quick identification of SNP shared by all replicates. (Variants summarization)
- Quick identification of Structural variations shared by all replicates. (Variants summarization)
- Summaries of all variation present inside a number of replicates. (Variants summarization)
- Obtain an indicator that summarize how similar two replicates are. (Files comparison)

This program seek .vcf files inside a folder (see `-p` option) and its sub folders. Files are grouped using theirs names or theirs parent folder (see `-s` option).
Each vcf of each group is compared with other vcf of the same group. The type of comparison depend on the option choose (see `-g` and `-b`). See “Variants summarization” and “ Files comparison” sections.

In order to run this program, download `main.sh`, `scan.py` and `compare.py` inside the same directory and use `bash main.sh -h` inside a linux terminal.

## Variants summarization:
### Purpose and functioning
The purpose of this mod is to summarize variations contained by a group of files. This is the default mod.

When this mod is used the program will seek each variant at each position of each file. Each position related to a variant shared by at least two replicates is displayed. Each displayed line contain a score that represent the proportion of file that contain a variant at this position and some other information (see below).

For each set of variant, a header mark by `###` is displayed. This header contain the name of this set of variant and various others information related to program parameters. Those parameters includes options used when the script was executed and files' versions.

Results are displayed on 6 Columns:
- SCORE : A percentage that represent the frequency of this variant. (Equal to GF / GM * 100)
- CHROM : The chromosome that contain this variant
- POS : The place of the variant on CHROM
- GF : A number that represent how many replicates contain a varaiant at this position.
- GM : The maximum number of times that this variant can be found inside this set of replicates. (equivalent to the number of replicates)
- OCUR : Store all variations with the following format : Variation=number of time that this variation have occured. Note that the sum off all occurrence from all variation can be greater than GM. This is due to the fact that multiple variation can be stored at a unique position in VCF files.


### Exploiting result
- The greater the percentage in the SCORE column, the more variant  at this position are shared between each replicates.
- Multiple items in the OCUR column mean that all variant a this position does not share the same variation. (This do reduce SCORE value.)
- The greater the OCUR is, the more alterations at this places can be found inside all replicates

Use the SCORE column to identity SNP and Structural variations shared by a significative number of replicates.

All lines inside this tab represent at least one variant that is share by at least two replicates.

## Files comparison :
### Purpose and functioning
The purpose of this mod is to make the identification of aberrant file easier.  

You can activate this mod with the option `-e` when you call `main.sh`.

When this mod is activated, the program will seek each combinations of files and look how many variants inside the first file match with variants inside the second file. Each file have two main score attributed to them:
- A percentage of similarity (GSCORE). This represent to what extent the variants inside this files end up in other files and to what extent variants of other files end up in this file.
- A score of inclusion (ISCORE). This score is also percentage. This percentage represent to what extent variants found inside this file can be found inside others files.

At the end of the program a the mean of all GSCORE files is calculated. Furthermore, since each combinations of files is compared a local ISCORE and a local GSORE is calculated for each files.

For each set of variant, a header mark by `###` is displayed. This header contain the name of this set of variant and various others information related to options used when the script was executed. This include execution options and about files' versions.

For each replicates a sub-header mark by `##` is displayed.  This sub-header contain the name of this replicate along with replicate’s global score (mean of all GSCORE) and with replicate’s inclusion score (mean of all ISCORE)

Under each replicates a tab is displayed:
- GSCORE :  Represent to what degree positions inside this replicate and FILE are similar. Equal to GF / GM * 100.
- GF : How many position inside this replicate have found at least one similar variant in FILE and vice versa.
- GM : How many position inside this replicates and inside FILE exist.
- ISCORE :  Represent to what degree positions inside this replicate are included inside FILE. How many position inside this replicate have found at least one similar variant in FILE. Equal to IF / IM * 100.
- IF  : How many position inside this replicate have found at least one similar variant inside FILE.
- IM : How many position inside this replicate exist.
- FILE : File name of the file that is compared to this replicate.

### Exploiting result
- Replicates with low GSCORE and high ISCORE might be replicates with a low number of positions. Those positions are likely to be contained in others replicates.
- Replicates with a low GSCORE and low ISCORE are really different from other replicates. They might be aberrant replicates
- Groups with low GSCORE contain poor quality replicates and / or aberrant replicates.
- Groups with high GSCORE  are likely to be composed of high quality replicates.
- Files with high ISCORE are likely to be a subset of other files.

# Main options
List of options accepted by `main.sh`.
- h) Return the help of `main.sh`.
- d) By default, files are opened during the search to verify that they are variant call format and not Vcard files. This also exclude multiple sample files (which are not supported). Use this option to turn it off.
- v) This program will display information during the process (file found, handled errors, progress bar...).
- g) This program will only return Files comparison. If let unspecified (and `-b` unspecified too), Variants summarization is returned.
- b) This program will return both Files comparison score and Variants summarization.
- c) Show files with their complete path.

- p) A path to folder. All files inside this folder, its sub-folders, its sub-sub-folder and so on will be passed in review. All .vcf files are used by this program. If let unspecified, ‘~’ is used.
- s) A Separator that will be used to group files (can not be 'none'). If unspecified parent folder will be used to group files. Here an example with “-” as a separator:
  - a file named P15-1.vcf will be inside the group “p15”,
  - a file named P30-1.vcf will be inside the group “p30”,
  - a file named P15-1-1.vcf will be inside the group “p15”,
  - a file named -P151.vcf will be inside the group “GroupNameLessFiles”,
  - a file named P15.vcf will be inside the group “SeparatorLessFiles”,
- o) How close two position should be to be compared to each-others. By default 0 is used. Positions that match together due to the offset gain half a similarity point. Offset can not cross chromosomes.
- t) Integer between 0 and 100. When two sequences are compared, they must have an alignment score greater or equal to this threshold to be considered identical. If threshold is unspecified, two sequences are considered similar if they are identical.
- r) A path toward a file. Result of these comparisons will be stored inside this file. If this file exist, it will be overwright. If unspecified result will be printed inside the console.

# How comparisons works
A position is the emplacement of a variant inside a genome.
Two position are considered similar when:
- Positions are close enough : positionB ∈ [positionA - offset ; positionA + offset]
- Positions have at least one alteration in common. When multiple alteration are present at the same position, all combinations between the two files are test and at least one combinations has to fulfill one of the following rules :
  - The two positions have a DEL alteration.
  - The two positions have a INS alteration.
  - The two positions have a DUP alteration.
  - The two positions have the exact same sequence (`-t` is unspecified).
  - The two positions have similar enough sequence (`-t` is an integer). We use a Smith-Waterman Algorithm.

# Programs improvement point and flaws
## Improvement point
- The process of grouping files could be more effective if file indexing and file grouping was done at the same time.
- The function `is_variant_call_format()` is costly in time when `-d` is let unspeciefied since it require to read each line of the files. Moreover, this function only seek for columns' legends and so is quite easy to fool.
- `scan.py` gain in usability if it could be called more easily from a linux terminal (using `getopt` library (this library was not allowed for this project)).
- It could a good idea to transform `compare_replicat()` into a generator in order to save memory.
- Add a sort option to sort result of Variants summarization.
- Score matrix for the Smith-Waterman Algorithm is hard coded.

## Known flaws
- VCF with multiple samples are not supported.
- Two position are considered equivalent when two DEL, two INS or two DUP are at the same position with no consideration for the length or for the sequence.
- Variant at the same position are stored together.

# Dependency
`python3`, `os` and `sys` library

# About this project
This project has been realized during the first semester of my master's degree in bio-informatics (initially I’m a biologist) at the university of Montpellier (France). The goal was to make a program to compare a number .vcf files. The only libraries authorized were `sys`, `os` and `re`. Custom objects (`class`) wasn’t authorized. 

Since this project was made during class to allow my teacher to evaluate me, I do not think I will maintain it unless I need it or unless someone ask me.
