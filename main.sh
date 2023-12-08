# ================================ HEADER ================================
# Requirements  : Python3, scan.py and compare.py
# Author        : Marchal Florent
# Copyright     : Copyright 2023, Marchal Florent
# License       : CC-BY-SA-4.0
# Credits       : Marchal Florent, Fiston-Lavier Anna-Sophie, Jacot Vetea
# Maintainer    : Marchal Florent
# Email         : flo.marchal2002@gmail.com
# Status        : Production
# ================================ ====== =================================

module_help="
This program will help you to:
    - examine whether or not replicates of an experiment are similar enough
    - extract a consensus of alternative sequences at set position.

This file use two python files (scan.py and compare.py) to seek, group and compare .vcf files.
A percentage is assigned to each group of vcf files. This score represent how similar files are inside this group.

/!\\ Does not support multi-samples files. /!\\

Options :
- h) Return the help of main.sh.
- d) By default, files are opened during the search to verify that they are variant call format and not Vcard files. This also exclude multiple sample files (which are not supported). Use this option to turn it off.
- v) This program will display information during the process (file found, handled errors, progress bar...).
- g) This program will only return Files comparison. If let unspecified (and b unspecified too), Variants summarization is returned.
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
- r) A path toward a file. Result of these comparisons will be stored inside this file. If unspecified result will be printed inside the console.
How this program work:
  1) Folders (-p) are scanned in order to find each variant call format files.
  2) Each variant call files are grouped using a separator or their parent folder (-s)
  3) Inside each groups, vcf files are compared to each others. When two files are compared, each positions of each files are passed in review.

Two position (chromosome and place on this chromosome) are considered similar when:
	- Positions are close enough : positionB ∈ [positionA - offset ; positionA + offset]
	- Positions have at least one alteration in common. When multiple alteration are present at the same position,
	  all combinations between the two files are test and at least one combinations has to fulfill one of the following rules :
    		- The two positions have a '<DEL>' alteration
    		- The two positions have a '<INS>' alteration
   	 	  - The two positions have a '<DUP>' alteration
    		- The two positions have the exact same sequence (<sequence_threshold> is None)
    		- The two positions have similar enough sequence (<sequence_threshold> is an integer)

--- --- --- Files comparison (g or b) --- --- ---
Use this mod when you want to known at which extend file inside a group are to similar together.
Two scores are given to each couple of files. The first score (GSCORE) represent how similar positions in the two files
are and the second one to what extent positions in the first file is included in the second.

For each set of variant, a header mark by ### is displayed. This header contain the name of this set of variant and various others information related to options used when the script was executed. This include execution options and about files' versions.

For each replicates a sub-header mark by ## is displayed.  This sub-header contain the name of this replicate along with replicate’s global score (mean of all GSCORE) and with replicate’s inclusion score (mean of all ISCORE)

under each replicates a tab is displayed:
  - GSCORE :  Represent to what degree positions inside this replicate and FILE are similar. Equal to GF / GM * 100.
  - GF : How many position inside this replicate have found at least one similar variant in FILE and vice versa.
  - GM : How many position inside this replicates and inside FILE exist.
  - ISCORE :  Represent to what degree positions inside this replicate are included inside FILE. How many position inside this replicate have found at least one similar variant in FILE. Equal to IF / IM * 100.
  - IF  : How many position inside this replicate have found at least one similar variant inside FILE.
  - IM : How many position inside this replicate exist.
  - FILE : File name of the file that is compared to this replicate.

--- --- --- Variants summarization (default or b)--- --- ---
Use this mod when you want to determine a consensus of variation inside a group .
Display a list of position sorted by similarity along with all ALT found at this place.

For each set of variant, a header mark by ### is displayed. This header contain the name of this set of variant and various others information related to program parameters. Those parameters includes options used when the script was executed and files' versions.

Results are displayed on 6 Columns:
  - SCORE : A percentage that represent the frequency of this variant. (Equal to GF / GM * 100)
  - CHROM : The chromosome that contain this variant
  - POS : The place of the variant on CHROM
  - GF : A number that represent how many time this variant is found inside this set of replicates. If a variant is found due to an offset, this score is raised by 0.5 instead of 1
  - GM : The maximum number of times that this variant can be found inside this set of replicates. (equivalent to the number of replicates)
  - OCUR : Store all variations with the following format : Variation=number of time that this variation have occured. Note that the sum off all occurrence from all variation can be greater than GM. This is due to the fact that multiple variation can be stored at a unique position in VCF files.
"
if [[ $# < 1]]
then
  echo "$module_help"; exit
fi

folder_path=~    # Path to the folder where scan.py will search vcf files.
separator=none   # A Separator to group vcf files. If empty, files wille be grouped using their parent folders.
offset=0         # Maximal distance between two position that can be compared to each others.
threshold=none   # threshold distance between two position that can be compared to each others.
open_files=true  # Do VCF files are opened in order to verify that they are not Vcard files.
quiet=true       # Do information about scan and comparison advancement are displayed.
output_file=none  # Do result are printed inside the console or write inside a file.
output_type=position
complete_names=false

while getopts 'hgbdcqp:s:o:t:r:' option;
do
  case "$option" in
  h) echo "$module_help"; exit
  ;;
  p) folder_path=$OPTARG
  ;;
  s) separator=$OPTARG
  ;;
  o) offset=$OPTARG
  ;;
  t) threshold=$OPTARG
  ;;
  d) open_files=false
  ;;
  v) quiet=false
  ;;
  r) output_file=$(readlink -e $OPTARG)
  ;;
  g) output_type="file"
  ;;
  b) output_type="both"
  ;;
  c) complete_names=true
  esac
done

folder_path=$(readlink -e $folder_path)
python3 scan.py $folder_path $separator $offset $threshold $open_files $quiet $output_file $output_type $complete_names
