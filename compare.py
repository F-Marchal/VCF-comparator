# encoding=utf-8
"""This file contain a number of functions that can be used to load and compare variant call format files.
/!\\ Does not support multi-samples files. /!\\

Critics :
    - A lot of improvement can be done on <_compare_position_alt> since two position are considered
        equivalent when two "<DEL>", two "<INS>" or two "<DUP>" are at the same position with no consideration for the
        length or for the sequence.
    - No functions to support multi-samples files
    - It could a good idea to transform <compare_replicat> into a generator in order to save memory
"""
__author__ = "Marchal Florent"
__copyright__ = "Copyright 2023, Marchal Florent"
__credits__ = ["Marchal Florent", " Fiston-Lavier Anna-Sophie", "Berard Severine"]
__license__ = "CC-BY-SA-4.0"
__version__ = "1.0.2"
__maintainer__ = "Marchal Florent"
__email__ = "flo.marchal2002@gmail.com"
__status__ = "Production"


def load_vcf_positions(path: str, keep_header: bool = False, keep_path: bool = False,
                       **line_options,
                       ) -> dict:
    """Load inside a dictionary all positions of a variant call format file. One position can be the key to multiple
    items as positions aren't uniq inside these kind of file.
    This function use <parse_vcf_line>.

    :param str path:            A path that lead to a Variant Call Format file
    :param bool keep_header:    Do the header is stored inside a list inside the dictionary.
                                    True: The dict will contain a key named 'header'. This key lead to a list.
                                    False: The header is not stored.
    :param bool keep_path:      Do <path> is stored inside the dictionary.
                                    True: The dict will contain a key named 'path'. This key lead to  <path>.
                                    False: The path is not stored.
    :param line_options:    keys for <parse_vcf_line>. Expected arguments :
        all_:       If True: None and undefined arguments are considered True
                    If False: None and undefined arguments are considered False
        parse_info: Does the column "INFO" is turned intoo a dictionary. (boolean)
        id_:        Do the "ID" column is included in the result ? (boolean)
        ref:        Do the "REF" column is included in the result ? (boolean)
        alt:        Do the "ALT" column is included in the result ? (boolean)
        qual:       Do the "QUAL" column is included in the result ? (boolean)
        filter_:    Do the "FILTER" column is included in the result ? (boolean)
        info:       Do the "INFO" column is included in the result ? (boolean)
        format_:    Do the "FORMAT" column is included in the result ? (boolean)
        samples:    Do all "SAMPLES" columns are included in the result ? (boolean) "SAMPLES" are stored inside a list.
    :return dict: {(chromosome, line_position (int)) : list of dict (dict from <parse_vcf_line>),
                   "header": list of lines inside the header,
                   "path": <path>}
    """
    # Preparation
    file = open(path)
    vcf_dictionary = {}

    # Prepare header and path if required.
    if keep_header:
        vcf_dictionary["header"] = []
    if keep_path:
        vcf_dictionary["path"] = path

    # The pos argument can not be passed to <line_options> by the user since this argument has to be true.
    if "pos" in line_options:
        raise TypeError("Unexpected keyword argument : 'pos'.")
    if "chrom" in line_options:
        raise TypeError("Unexpected keyword argument : 'chrom'.")

    # Start file loading.
    for i, lines in enumerate(file):
        # Remove line break
        while lines[-1] == "\n":
            lines = lines[:-1]

        if keep_header and lines[0:2] == "##":
            # Header lines are characterised by "##"
            vcf_dictionary["header"].append(lines)

        elif lines[0:1] == "#":
            # Legend lines is characterised by "#"
            continue

        else:
            # Extract line
            try:
                line = parse_vcf_line(lines, pos=True, chrom=True, **line_options)
            except IndexError as E:
                raise IndexError(f"{E} Line {i}")

            # Convert POS
            try:
                position = int(line["POS"])
            except ValueError as E:
                raise ValueError(f"Can not turn position into an integer (line {i}, file {path}). Line ignored")

            # If we forget CHROM, two position on two chromosome can be treated as if there were at the same place.
            position = (line["CHROM"], position)

            # Save the results
            if position not in vcf_dictionary:
                vcf_dictionary[position] = [line]
            else:
                vcf_dictionary[position].append(line)

    file.close()
    return vcf_dictionary


def parse_vcf_line(line: str, parse_info: bool = True, all_=True,
                   chrom=None, pos=None, id_=None, ref=None, alt=None, qual=None, filter_=None, info=None,
                   format_=None, samples=None,
                   ) -> dict:
    """Take a string that represent a body line of a variant call format file and turn it into a dictionary.

    :param str line:        A line from the body of a variant call format file
    :param parse_info:      Does the column "INFO" is turned intoo a dictionary.
    :param bool all_:       If True: None arguments are considered True
                            If False: None arguments are considered False
    :param bool chrom:      Do the column "CHROM" is included in the result ?
    :param bool pos:        Do the column "POS" is included in the result ?
    :param bool id_:        Do the column "ID" is included in the result ?
    :param bool ref:        Do the column "REF" is included in the result ?
    :param bool alt:        Do the column "ALT" is included in the result ?
    :param bool qual:       Do the column "QUAL" is included in the result ?
    :param bool filter_:    Do the column "FILTER" is included in the result ?
    :param bool info:       Do the column "INFO" is included in the result ?
    :param bool format_:    Do the column "FORMAT" is included in the result ?
    :param bool samples:    Do all "SAMPLES" columns are included in the result ? "SAMPLES" are stored inside a list.
    :return:
    """
    # Remove line breaks for the line and split the line in columns
    while line[-1] == "\n":
        line = line[:-1]
    split_line = line.split("\t")

    # Lines inside .vfc files are supposed to have at least 8 columns.
    if len(split_line) < 8:
        print(split_line, line.replace("\t", ","))
        raise IndexError(f"Can not extract 8 fields from '{line}'. All fields has to be splited by \\t (tabulations).")

    # Fill <line_dict> with the correct columns.
    line_dict = {}
    if _test_in_PCL(all_, chrom): line_dict["CHROM"] = split_line[0]
    if _test_in_PCL(all_, pos): line_dict["POS"] = split_line[1]
    if _test_in_PCL(all_, id_): line_dict["ID"] = split_line[2]
    if _test_in_PCL(all_, ref): line_dict["REF"] = split_line[3]
    if _test_in_PCL(all_, alt): line_dict["ALT"] = split_line[4]
    if _test_in_PCL(all_, qual): line_dict["QUAL"] = split_line[5]
    if _test_in_PCL(all_, filter_): line_dict["FILTER"] = split_line[6]
    if _test_in_PCL(all_, info): line_dict["INFO"] = parse_vcf_line_info(split_line[7]) if (parse_info
                                                                                            is True) else split_line[7]
    if _test_in_PCL(all_, format_): line_dict["FORMAT"] = split_line[8] if len(split_line) > 8 else None
    if _test_in_PCL(all_, samples): line_dict["samples"] = split_line[9:] if len(split_line) > 9 else []

    return line_dict


def _test_in_PCL(all_: bool, item: bool = None) -> bool:
    """Internal function. Function used by <parse_vcf_line> to determine if an information will be inside
    the returned dict.
    :param boolean all_:   If True: None <item> is considered True
                           If False: None <item> is considered False
    :param boolean item: True, False or None
    :return bool: Follow the following tab:
            item    all     result
            None    True    True
            None    False   False
            False   True    False
            False   False   False
            True    True    True
            True    False   True
    """
    return item is True or (all_ is True and item is not False)


def parse_vcf_line_info(info: str) -> dict:
    """Store "INFO" column inside a dictionary all info are used as key and theirs attached values are used as value.
    :param info:    String from the INFO column of a line.
    """
    intel = {}
    # All "INFO" are delimited by a ";".
    for items in info.split(";"):
        if "=" not in items:
            # Save this item
            intel[items] = ""
            continue

        # Extract key and value and save them inside <intel>
        key, value = items.split("=", 1)
        intel[key] = value

    return intel


def compare_replicat(offset: int = 0, sequence_threshold: float = None, quiet: bool = True,
                     **replicates) -> (dict[str], dict[tuple[int, str]]):
    """Compare a number of replicates using their positions alterations. All replicates are compared two per two.
    A score of global similarity and a score of inclusion is given for all replicates.
    Also, a summary of which position are the most common is returned.

    --- --- Replicates comparisons --- ---
        Let 'a' a replicate, 'b' a second replicate.
        The Score of global similarity is calculated with the following formula :
            ((number of position 'a' that match with position inside 'b') + (number of position 'a' that match with
            position inside 'b')) / ((number of position in 'a') + (number of position in 'b'))

        The Score of inclusion is calculated with the following formula :
            (number of position 'a' that match with position inside 'b') / (number of position in 'a')

    --- --- Positions comparisons --- ---
        Let 'a' a replicate, 'b' a second replicate.
        When a position inside 'a' match with a position inside 'b' those positions are stored inside a dict that
        summarize how many times those position has a match and which alternative has been found.
        Alternative are store inside a dictionary and are used as key. Those key lead to a counter of occurrences.
        When two position match due to an offset, this counter is increased by 0.5 and when they match without any
        offset this counter is increased by 1. Since multiple alterations can touch a position at the same time, this
        counter isn't limited by the number of replicates.

    --- --- Rules for a match --- ---
    Two position can match together if :
        - Positions are close enough : positionB ∈ [positionA - offset ; positionA + offset]
        - Positions have at least one alteration in common:
            - The two positions have a "<DEL>" alteration
            - The two positions have a "<INS>" alteration
            - The two positions have a "<DUP>" alteration
            - The two positions have the exact same sequence (<sequence_threshold> is None)
            - The two positions have similar enough sequence (<sequence_threshold> is an integer)

    -- --- Functions arguments --- ---
    :param int offset:          How close the position should be to be compared to each-others.
    :param int sequence_threshold:  How similar two sequences must be to consider them  as identical.
                                    If None, sequences has to be the same.
    :param bool quiet:          False: This function will print a progress bar to show the progression.
                                True: This function will not print anything.
    :param replicates: At least two replicates. replicates names can not be '__MEANS__'.
        replicate_name=replicate_dict.
    :return tuple[dict]:
        - score_dict = {
            "__MEANS__": {"__MEANS__": [Global mean],
                          "<replicate_name>": [Global mean, inclusion mean]
            },
            "<replicate_name>": {
                                "<other_replicate_name>" : (global score,
                                                            number of global match,
                                                            maximum match number,
                                                            inclusion score (of <replicate_name> in <other_replicate_name>),
                                                            number of inclusion match,
                                                            maximum inclusion
                                                            match number)
            },
        }
        - position_dict = {
                (chromosome (str), position (int)): [number of replicates, {
                                    [ALT at these positions: occurrences (int)]
                                }],
        }
    """
    # Some verification
    if len(replicates) < 2:
        raise ValueError("Not enough replicate provided. At least two replicate are expected.")
    if "__MEANS__" in replicates:
        raise NameError("Can not compute replicates with '__MEANS__' as name.")
    if (sequence_threshold is not None) and (not 0 <= sequence_threshold <= 100):
        raise ValueError("<sequence_threshold> should be None or a number "
                         "greater or equal to 0 and lower or equal to 100")
    if offset < 0:
        offset = 0

    # Prepare some variable for __MEANS__ and for the progressbar
    number_of_replicates = len(replicates)
    number_of_comparison_per_replicates = number_of_replicates - 1
    number_of_comparison = 0
    for i in range(0, number_of_replicates):
        number_of_comparison += i

    # Prepare the progress bar
    progress_bar = 0            # Items displayed inside the progress bar.
    progress_per_replicate = 1  # Amount of progression granted when two replicates are compared.
    total_progress = 0          # Total of progression.
    if not quiet:
        progress_per_replicate = len("|----------|----------|----------|----------|") / number_of_comparison
        print("           0          25         50         75        100")
        print("Progress : |----------|----------|----------|----------|")
        print("           ", end="", flush=True)

    # Prepare the final_dictionary and replicates:
    score_dict = {"__MEANS__": {"__MEANS__": [0]}}
    positions_dict = {}
    replicates_list = [(dict_name, dic_, len(dic_)) for dict_name, dic_ in replicates.items()]
    comparison_errors = []

    # Begin the comparisons:
    for i, (main_name, main_dict, main_length) in enumerate(replicates_list[:-1]):
        for second_name, second_dict, second_length in replicates_list[i + 1:]:
            # -- -- Loop Breakdown -- --
            # Let 'a' a replicate, 'b' a second replicate, 'max' the number of items inside <replicates_list>,
            # 'n' a number between 0 and 'max' - 1.
            #
            # This two nested loops work by comparing the nth item with all items between 'n' + 1 and 'max'.
            # Since a comparison of 'a' with 'b' is equal to a comparison of 'b' with 'a', this loops can fill
            # the result of 'a' with 'b' and the result of 'b' with 'a' at the same time.
            # Effectively, the last item of <replicates_list> has been compared with all other replicates when the first
            # loop reach his position.
            # The complexity of this nested loops follow a binomial coefficient of parameter ('max', 2)
            # -- -- -- -- -- -- -- -- --
            # <main_match> and <second_match> will respectively contain all match made by main_dict's and second_dict's
            #           Positions. We use a set instead of a simpler counter in order to avoid counting
            #           multiple times one match when <offset> is greater than 0.
            # -- -- -- -- -- -- -- -- --

            # Comparisons variable
            main_match = set()
            second_match = set()
            max_matches = main_length + second_length

            # Match finder
            for initial_pos in main_dict:  # <initial_pos> is a position without any offset
                for j in range(-offset, offset + 1):
                    current_pos = (initial_pos[0], initial_pos[1] + j)
                    # <current_pos> is a position with an offset of <j>

                    if current_pos not in second_dict:
                        # No position match with this offset.
                        continue
                    else:
                        try:
                            results = _compare_position_alt(main_dict[initial_pos], second_dict[current_pos],
                                                            sequence_threshold=sequence_threshold)
                        except KeyError as E:
                            if not quiet:
                                comparison_errors.append(f"Can not proceed to the comparison of the position "
                                                         f"{initial_pos} (from {main_name}) "
                                                         f"with the position {current_pos} (from {second_name}) : {E}")
                            continue

                        if not results:
                            # This position does not fulfill requirement and so can not be kept as identical.
                            continue

                    # Store those positions as they passed earlier verification.
                    main_match.add(initial_pos)
                    second_match.add(current_pos)

                    # Store Positions
                    if initial_pos not in positions_dict:
                        initials_alterations = {}
                        positions_dict[initial_pos] = [set(), number_of_replicates, initials_alterations]
                    else:
                        initials_alterations = positions_dict[initial_pos][-1]

                    if current_pos not in positions_dict:
                        current_alterations = {}
                        positions_dict[current_pos] = [set(), number_of_replicates, current_alterations]
                    else:
                        current_alterations = positions_dict[current_pos][-1]

                    # <positions_dict[initial_pos][0]> is a set for the same raisons as <main_match>.
                    positions_dict[initial_pos][0].add(second_name)
                    positions_dict[current_pos][0].add(main_name)

                    for alt in results:
                        if isinstance(alt, tuple):
                            alt, second_alt = alt
                        else:
                            second_alt = alt

                        if alt not in initials_alterations:
                            initials_alterations[alt] = 0
                        if second_alt not in current_alterations:
                            current_alterations[second_alt] = 0

                        initials_alterations[alt] += 0.5
                        current_alterations[second_alt] += 0.5

            # Assure that names are in <score_dict>:
            if main_name not in score_dict:
                score_dict[main_name] = {}
                score_dict["__MEANS__"][main_name] = [0, 0]
            if second_name not in score_dict:
                score_dict[second_name] = {}
                score_dict["__MEANS__"][second_name] = [0, 0]

            # Global comparison score
            global_match = len(main_match) + len(second_match)
            if max_matches != 0:
                global_percent = round(global_match / max_matches * 100, 2)
                global_result = (global_percent, global_match, max_matches)
            else:
                global_result = (100, 0, 0)

            # Save results in <score_dict>:
            second_score = (
                round(len(main_match) / main_length * 100, 2) if main_length > 0 else 100,
                len(main_match),
                main_length)

            main_score = (
                round(len(second_match) / second_length * 100, 2) if second_length > 0 else 100,
                len(second_match),
                second_length)

            score_dict[main_name][second_name] = (*global_result, *second_score)
            score_dict[second_name][main_name] = (*global_result, *main_score)

            # Save results in __MEANS__
            score_dict["__MEANS__"][main_name][0] += global_result[0] / number_of_comparison_per_replicates
            score_dict["__MEANS__"][main_name][1] += main_score[0] / number_of_comparison_per_replicates
            score_dict["__MEANS__"][second_name][0] += global_result[0] / number_of_comparison_per_replicates
            score_dict["__MEANS__"][second_name][1] += second_score[0] / number_of_comparison_per_replicates
            score_dict["__MEANS__"]["__MEANS__"][0] += global_result[0] / number_of_comparison

            # Update the progress bar
            if not quiet:
                total_progress += progress_per_replicate
                old_progress = progress_bar
                progress_bar = int(round(total_progress, 0))
                if (number_of_new_equals := progress_bar - old_progress) != 0:
                    print("=" * number_of_new_equals, end="", flush=True)

    if not quiet:
        # End the progress bar
        print()
        for items in comparison_errors:
            # Show problematics lines
            print(items)
        print()

    # End the function
    return score_dict, positions_dict


def _compare_position_alt(main_items: list, second_items: list, sequence_threshold: float = None):
    """Compare two list of alteration to say if at least one item of <main_items> match with at least one item inside
    <second_items>.

    A match is count if :
    - At least one alteration inside <main_items> and at least one alteration inside <second_items> is a "<DEL>". (1)
    - At least one alteration inside <main_items> and at least one alteration inside <second_items> is a "<INS>" (1)
    - At least one alteration inside <main_items> and at least one alteration inside <second_items> is a "<DUP>" (1)
    - At least one alteration inside <main_items> and at least one alteration <second_items> are an identical
        sequence (<sequence_threshold> is None)  (2)
        OR are identical enough (percent_alignment >= <sequence_threshold>). (3)
    - The two positions have similar enough sequence (<sequence_threshold> is an integer)

    :param list[dict] main_items:      A list of lines from a vcf file (dictionaries with "ALT" inside)
    :param list[dict] second_items:    Another list of lines from a vcf file (dictionaries with "ALT" inside)
    :param float sequence_threshold:  How similar two sequences must be to consider them identical.
                                      If None, sequences has to be the same.
    :return list: a list of all match that as occurred between <main_items> and <second_items>.
    """
    matches = []
    i = 0
    while i < len(main_items):
        # Update variables
        main_variations = main_items[i]
        alt = main_variations["ALT"]
        j = 0

        while j < len(second_items):
            # Update variables
            second_variations = second_items[j]
            second_alt = second_variations["ALT"]

            # Comparison
            if alt in ("<DEL>", "<INS>", "<DUP>", "H") or second_alt in ("<DEL>", "<INS>", "<DUP>", "H"): # (1)
                if alt == second_alt:
                    matches.append(alt)

            elif sequence_threshold is None:  # (2)
                if alt == second_alt:
                    matches.append(alt)

            else:   # (3)
                percentage = seq_percent_alignment(seqA=alt, seqB=second_alt)
                if percentage >= sequence_threshold:
                    matches.append((alt, second_alt))

            j += 1
        i += 1
    return matches


def seq_percent_alignment(seqA: str, seqB: str, gap: int = 3, substitution: dict = None, max_score=None) -> float:
    """Smith-Waterman Algorithme that compare two sequence and return a percent of similarity.
    Gap score and scores inside the substitution matrix have to be positive or null. Otherwise, this function can
    not return a correct percentage.

    :param str seqA:            A sequence to analise
    :param str seqB:            A second sequence to analise
    :param int gap:             Score of a gap (must be positive or null)
    :param dict substitution:   Dict that contain all possible substitutions and give them a score. This score must be
                                positive or null
                                If None, the following matrix is used:
                                       A   C   G   T
                                    A  9   0   4   0
                                    C  0   9   0   4
                                    G  4   0   9   0
                                    T  0   4   0   9
    :param int max_score:       The maximum score that can be found in <substitution>.
                                    If None : <substitution>['A']['A'] is used
    :return float: percent of similarity. (greater or equal to 0 and lower or equal to 100)
    """
    number_of_lines = len(seqA) + 1
    number_of_columns = len(seqB) + 1

    if substitution is None:
        # Generate a substitution matrix
        substitution = {
            "A": {
                "A": 9,
                "C": 0,
                "G": 4,
                "T": 0,
            },
            "C": {
                "A": 0,
                "C": 9,
                "G": 0,
                "T": 4,
            },
            "G": {
                "A": 4,
                "C": 0,
                "G": 9,
                "T": 0,
            },
            "T": {
                "A": 0,
                "C": 4,
                "G": 0,
                "T": 9,
            },
        }
    if not max_score:
        # Find a value for max_score
        max_score = substitution["A"]["A"] * max(number_of_lines - 1, number_of_columns - 1)

    # Matrix generation
    empty_line = [0 for j in range(0, number_of_columns)]
    matrix = [empty_line[:] for i in range(0, number_of_lines)]

    # Initialisation
    for i in range(0, number_of_lines):
        matrix[i][0] = i * gap
    for j in range(0, number_of_columns):
        matrix[0][j] = j * gap

    # Fill the matrix
    for i in range(1, number_of_lines):
        for j in range(1, number_of_columns):
            try:
                score_sub = matrix[i-1][j-1] + substitution[seqA[i-1]][seqB[j-1]]
            except KeyError as E:
                raise KeyError(f"Item not found inside the substitution matrix : {E}")

            score_insert = matrix[i][j-1] + gap
            score_del = matrix[i-1][j] + gap
            matrix[i][j] = max(score_sub, score_insert, score_del)

    # Return the score
    return matrix[-1][-1] / max_score * 100

