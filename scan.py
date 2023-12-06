# encoding=utf-8
"""This file contain a number of functions that can be used to 	identify and group .vcf files. Variant call format files
can be grouped using their folder <group_file_by_folder> or using their names and a separator <group_file_by_name>.
This file also contain a <main> that index and compare every .vcf files present below a path.

/!\\ Does not support multi-samples files. /!\\

If you decide to call this file from Bash, here a list of accepted arguments (see <main>):
    1 - folder_path
    2 - offset
    3 - threshold
    4 - open_files
    5 - quiet
    6 - output_file
    7 - complete_names
    8 - output_type

Critics:
    - The process of grouping files could be more effective if file indexing and file grouping was
        done at the same time.
    - The function <is_variant_call_format> is costly in time when open_file=True since it require to read each line of
        the files. Moreover, this function only seek for columns' legends and so is quite easy to fool.
    - This file could gain in usability if it could be called more easily from Bash (using getopt) (this library was
        not allowed for this project).
    - A time counter would have been a good idea inside <find_variant_call_format_file> to show to the user that
        the program isn't locked. (Indexing files can be quite long.) (time library was not allowed for this project)
"""

__author__ = "Marchal Florent"
__copyright__ = "Copyright 2023, Marchal Florent"
__credits__ = ["Marchal Florent", " Fiston-Lavier Anna-Sophie"]
__license__ = "CC-BY-SA-4.0"
__version__ = "1.0.1"
__maintainer__ = "Marchal Florent"
__email__ = "flo.marchal2002@gmail.com"
__status__ = "Production"


import os
import sys
import compare


def is_variant_call_format(path: str, open_file: bool = True, one_sample_only: bool = True, quiet: bool = True) -> bool:
    """This function tell if a file (determined by a path (<path>)) is a variant cell format.
    This function need to open (<open_file>=True) the file in order to discriminate variant cell format and vCard and
    in order to discriminate files with multiple samples.

    :param str path:    A path to find to the targeted file.
    :param open_file:   Do we look what is inside this file ? default=True (discriminate variant cell format and vCard)
    :param one_sample_only:    If True, only file with one sample inside are used. REQUIRE <open_file>=True.
    :param bool quiet:  If False .VCF that does not respect requirement are showed.
    :return bool:       Do this file is a variant_call_format ?"""
    extension = path.split(".")[-1]
    if extension.lower() != "vcf":
        return False
    elif open_file is not True:
        # .vcf file.
        return True

    # Define a marker in order to determine whether this file is a variant call format file.
    file_marker = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    file_marker_length = len(file_marker)

    # <open_file> is True and file is .vcf. Let's open it.
    result = False
    try:
        file = open(path)
    except PermissionError as E:
        if not quiet: print(f"Can not open {path} ({E})")
        return result

    # Loop through this file until <file_marker> is found.
    line = "temp"
    while result is False and line:
        result = line[0:file_marker_length] == file_marker

        if not result:
            line = file.readline()
        # Else : The loop will end soon (<result> is True) and we need to keep the line untouched in order to
        # verifications related to <one_sample_only>.

    if result is False:     # <file_marker> is not found.
        if not quiet: print(f"This file has .vcf extension but does not match with .vcf signature : {path}")

    elif one_sample_only is True and result is True and len(line.split("\t")) > 10:
        # If only one sample is inside this file then the legend can not be greater than 10 (8 required columns
        # + FORMAT + SAMPLE1)
        if not quiet: print(f"This file can not be used since it contain multiple samples : {path}")
        result = False

    file.close()
    return result


def find_variant_call_format_file(path: str, open_file: bool = True, quiet: bool = True,
                                  one_sample_only: bool = True) -> list:
    """List all vcf file inside a folder (<path>), its sub-folders, its sub-sub-folders and so on.
    This function use <is_variant_call_format> and so is able to discriminate variant call format from vCard.

    :param str path:           Path to a folder.
    :param bool open_file:     Do we open file to verify if files are variant call format ? (exclude vCard file).
    :param bool quiet:         If false, errors handle by this function are displayed.
    :param one_sample_only:    If True, only file with one sample inside are returned. REQUIRE <open_file>=True.
    :return list: A list of path that led to .vcf file that respect function's settings.
    """
    list_of_vcf = []
    if not path[-1] in ("/", "\\"):
        # Assure that <path> will be considered as a folder.
        path += "/"

    for items in os.listdir(path):
        item_path = path + items

        if os.path.isdir(item_path):
            try:
                list_of_vcf.extend(find_variant_call_format_file(item_path, open_file=open_file, quiet=quiet,
                                   one_sample_only=one_sample_only))
            except PermissionError as E:    # Might happen when a <path> next close to the root is given.
                if not quiet: print(f"Can not access to '{item_path}' ({E}).")

        elif is_variant_call_format(item_path, open_file, quiet=quiet, one_sample_only=one_sample_only):
            if not quiet: print("File found : " + item_path)
            list_of_vcf.append(item_path)

    return list_of_vcf


def group_file_by_name(list_of_file: list[str], separator: str = "-") -> dict[str:list[str]]:
    """Group files inside a dictionary of list using their names.
    Eg : list_of_file = ["P15-1", "P30-1", "P15-Alpha"] separator="-":
        -> {"P15": ["P15-1", "P15-Alpha"], "P30": ["P30-1"]}

    :param  list[str] list_of_file: A list of file path / names.
    :param str separator:           A separator that should be inside each file name
                                        If the separator isn't inside a name, files will be grouped
                                            as "SeparatorLessFiles"
                                        If there isn't any characters before the separator, files will be grouped
                                            as "GroupNameLessFiles"
                                        If separator is inside a name, files will be grouped
                                            as f"{String before <separator>}"
    :return dict:   A dictionary with the following form : {"{GroupName}" : ["{FilePath1}", "{FilePath1}", ...]}
    """
    groups = {}
    for file_path in list_of_file:
        # Remove path
        file_name = file_path.split("/")[-1]
        file_name = file_name.split("\\")[-1]

        # remove extension
        extension_less_name = ".".join(file_name.split(".")[:-1])  # Remove the extension but keep everything else

        split_name = extension_less_name.split(separator)
        if len(split_name) >= 2:
            # The splitter has been found, let's try to extract a group name.
            group_name = split_name[0]
            if len(group_name) == 0:
                # No group has been found (<file_name> start by <separator>)
                group_name = "GroupNameLessFiles"
        else:
            # Separator isn't inside <file_name>
            group_name = "SeparatorLessFiles"

        # Memorize this file
        if group_name in groups:
            groups[group_name].append(file_path)

        else:   # A new group is required
            groups[group_name] = [file_path]

    return groups


def group_file_by_folder(list_of_path: list[str]) -> dict[str:list[str]]:
    """Group files inside a dictionary using their path. File inside the same directory are grouped together.

    :param list list_of_path: A list of path.
    :return dict: A dictionary of list : {"{GroupName}" : ["{FilePath1}", "{FilePath1}", ...]}
    """
    groups = {}

    for paths in list_of_path:
        # <paths> cleaning
        if "\\" in paths:
            paths = paths.replace("\\", "/")

        # Extract a path toward the parent folder.
        if "/" in paths:
            split = paths.split("/")
            folder_path = "/".join(split[:-1])

        else:
            folder_path = "UNKNOWN"

        # Store files using their groups
        if folder_path in groups:
            groups[folder_path].append(paths)

        else:
            groups[folder_path] = [paths]

    return groups


def main(path, separator: str = "", offset: int = 0, threshold: float = None, open_files: bool = True,
         quiet: bool = True, output_file: str = None,
         complete_names: bool = False, output_type: str="position"):
    """Seek .vcf files inside a folder and its sub folders. Files are groups using theirs names or theirs parent folder.
    Each vcf of each group is compared with other vcf of the same group. A score of similarity is then displayed.
    See <compare.compare_replicat> to know how the score of similarity is determined.

    :param str path:            A path from where .vcf files are gathered.
    :param str separator:       A string that is used to split files names in order to groups
                                    them (see <group_file_by_name>). If this string is empty, files are groups
                                    using their parent folder (<group_file_by_folder>).
    :param int offset:          Define a range around all positions when two files are compared.
                                    Two position can be compared if the second one is inside the range of the first one.
    :param float threshold:     When two sequences are compared, they must have a percent of identity greater or
                                    equal to this threshold to be considered as identical.
                                    If threshold is None, two sequences are identical if they have the
                                        exact same sequence.
    :param bool open_files:     Does files are opened during the indexing period to identify if files are variant format
                                    call ? Use this if you have Vcard files in your computer
    :param bool quiet:          Do This program will print information during file processing ?
    :param str output_file:     A path toward a file where results can be saved. If None, results are printed
                                    intoo the consol.
    :param bool complete_names:       Do files' names are displayed using their full path
    :param str output_type:     'position', 'both' or 'file'. Define the type of data displayed.
                                    'position': Variants summarization
                                    'file': Files comparison
                                    'both': Both Files comparison and Variants summarization

    """
    # Make some verification
    if not os.path.isdir(path):
        raise NameError(f"Directory expected for <path>. Got : {path}")
    elif output_file is not None:
        if os.path.isdir(output_file):
            raise NameError(f"File expected for <output_file>. Got : {output_file}")
        output_file += ".txt"
    if output_type not in ("file", "both", "position"):
        raise ValueError(f"<output_type> is expected to be 'file', 'both' or 'position'. Got : {output_type}")

    # prepare some variables
    str_settings = f"path='{path}';separator='{separator}';offset={offset};threshold={threshold};"
    str_settings += f"open_files={open_files};quiet={quiet};output_file={output_file};complete_names={complete_names};"
    str_settings += f"output_type={output_type};SVersion={__version__};CVersion={compare.__version__}\n"

    if not quiet: print("settings : ", str_settings)

    # --- --- Find all vcf files --- ---
    if not quiet: print("======== Indexing =========\nIndexing .vcf files. This can take some time.")
    list_of_path = find_variant_call_format_file(path, open_files, quiet=quiet, one_sample_only=True)
    if not quiet: print(f"Indexing done : {len(list_of_path)} files found.")

    #  --- --- Group files --- ---
    if not quiet:
        p_sep = separator if len(separator) > 0 else "folder names"
        print(f"======== Group files =========\n"
              f"Separator used : {p_sep}")

    if len(separator) == 0:
        grouped_files = group_file_by_folder(list_of_path)
    else:
        grouped_files = group_file_by_name(list_of_path, separator=separator)

    if not quiet:
        print(f"{len(grouped_files)} Groups made :")
        for group_names, content in grouped_files.items():
            print(f"    {group_names} : {len(content)} item(s)")

    #  --- --- file Processing --- ---
    file_legend = "#GSCORE\tGF\tGM\tISCORE\tIF\tIM\tFILE\n"
    position_legend = "#SCORE\tCHROM\tPOS\tGF\tGM\tOCUR\n"
    for groups_name, list_of_files in grouped_files.items():
        if not quiet: print(f"===== Group : '{groups_name}' ====")

        # Groups can not be too smalls
        if len(list_of_files) < 2:
            if not quiet: print(f"'{groups_name}' group is too small : {len(list_of_files)} item(s) / {2}.")
            continue

        # --- Load files ---
        group_dict = {}
        for paths in list_of_files:
            # Handle errors raised by <compare.load_vcf_positions>
            try:
                vcf_dict = compare.load_vcf_positions(paths, all_=False, alt=True)
            except IndexError as E:
                if not quiet: print(f"Can not load {paths} : {E}")
                continue
            except ValueError as E:
                if not quiet: print(E)
                continue
            else:
                # Save results
                group_dict[paths] = vcf_dict

        # Groups can not be too smalls (again). Group can reduce in volume if some file can not be load.
        if len(group_dict) < 2:
            if not quiet: print(f"Can not use this group. Not enough file can be loaded.")
            continue

        # --- File comparisons ---
        score_dict, position_dict = compare.compare_replicat(offset=offset, sequence_threshold=threshold,
                                                             quiet=quiet, **group_dict)

        # --- Display results ---
        # Group header
        group_header_has_been_displayed = False
        paragraph = ""
        # --- Output scores related to files ---
        if output_type in ("file", "both"):
            paragraph = (f"###{groups_name}\tglobal={round(score_dict['__MEANS__']['__MEANS__'][0], 4)}%\t"
                         f"settings: {str_settings}")
            group_header_has_been_displayed = True

            for paths, comparisons in score_dict.items():
                if paths == "__MEANS__":
                    # This key is not a file. This key is used to store file's means.
                    continue

                # Reduce name size
                if complete_names is False:
                    name = paths.split("/")[-1].split("\\")[-1]
                else:
                    name = paths

                # Generate a text related to this file.
                path_means = score_dict['__MEANS__'][paths]
                paragraph += f"##{name}\tglobal={round(path_means[0], 4)}\tinclusion={round(path_means[1], 4)}\n"
                paragraph += file_legend

                # Assure that result are sorted by best score.
                sorted_comparisons = sorted(comparisons.items(), key=lambda item: item[1][0], reverse=True)

                # Display information related to each comparison
                for second_path, results in sorted_comparisons:
                    paragraph += "\t".join([str(items) for items in results]) + "\t" + second_path + "\n"

        # --- Output scores related to positions ---
        if output_type in ("position", "both"):
            if group_header_has_been_displayed is False:
                paragraph += (f"###{groups_name}\tsettings: {str_settings}")
                group_header_has_been_displayed = True
            paragraph += position_legend

            # Assure that positions are displayed from the greater occurrence to the lowest.
            sorted_positions = sorted(position_dict.items(), key=lambda item: round(item[1][0] / item[1][1] * 100, 4),
                                      reverse=True)

            # Display positions
            for position, (found, max_, elements) in sorted_positions:
                elements_detail = ";".join([f"{key}={item}" for key, item in elements.items()])
                chrom, pos = position
                paragraph += f"{round(found / max_ * 100, 4)}\t{chrom}\t{pos}\t{found}\t{max_}\t{elements_detail}\n"

        # Output result
        if output_file:
            with open(output_file, mode="a") as file:
                file.write(paragraph)
        else:
            print(paragraph)


if __name__ == "__main__":  # If this file isn't an import.
    sys_args = sys.argv[1:]
    args_length = len(sys_args)

    # Translate arguments from sys.argv[1:] and call <main>()

    # path
    if args_length >= 1 and sys_args[0] != "none":
        main_path = sys_args[0]
    else:
        main_path = os.path.expanduser("~")

    # separator
    if args_length >= 2 and sys_args[1] != "none":
        main_separator = sys_args[1]
    else:
        main_separator = ""

    # offset
    if args_length >= 3:
        try:
            main_offset = int(sys_args[2])
        except ValueError:
            raise ValueError(f"Integer expected for the 'offset' option. Got : {sys_args[2]}")
    else:
        main_offset = 0

    # threshold
    if (args_length >= 4) and sys_args[3] != "none":
        try:
            main_threshold = float(sys_args[3])
        except ValueError:
            raise ValueError(f"Float expected for the 'threshold' option. Got : {sys_args[3]}")
    else:
        main_threshold = None

    # open_files
    if args_length >= 5:
        if sys_args[4] in ("true", "1", "y"):
            main_open_files = True
        else:
            main_open_files = False
    else:
        main_open_files = True

    # quiet
    if args_length >= 6 and sys_args[5] in ("true", "1", "y"):
        main_quiet = True
    else:
        main_quiet = False

    # output file
    if args_length >= 7 and sys_args[6] != "none":
        main_output = sys_args[6]
    else:
        main_output = None
    
    # Output format
    if args_length >= 8 and sys_args[7] in ("file", "both", "position"):
        main_output_type = sys_args[7]
    else:
        main_output_type = "position"
        main_output_type = "file"
        
    # complete file name
    if args_length >= 9:
        if sys_args[8] in ("true", "1", "y"):
            main_complete_names = True
        else:
            main_complete_names = False
    else:
        main_complete_names = False
        
    # main
    main(
        path=main_path,
        separator=main_separator,
        offset=main_offset,
        threshold=main_threshold,
        open_files=main_open_files, 
        quiet=main_quiet, 
        output_file=main_output,
        output_type=main_output_type,
        complete_names=main_complete_names,
    )
