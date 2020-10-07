"""Takes a manifest file, and does all to all comparision to detect mosaic mutations"""
__author__ = "Vivekananda Sarangi"
__email__ = "sarangi.vivekananda@mayo.edu"
__status__ = "Development"

import os
import argparse
import sys
import Util
import gzip
import pandas as pd
import numpy as np
from math import sqrt
import pickle
import matplotlib.pyplot as plt
from subprocess import PIPE, Popen
from matplotlib.patches import Rectangle

class ALL2():
    def __init__(self):
        parser = argparse.ArgumentParser(description='All to all comparision',
                                         usage=''' python ALL2.py <command> [<args>]
        Three commands to use for SNVs and small INDELS:
                score --> Generates mosaic and germline scores.
                call --> Based on score cut-off generates sample level files/plots for mosaic,
                        germline mutations and plots variant allele frequency, mutation spectrum
                matrix --> Plot the mutation matrix
        Three commands to use for structural variants :
                score_sv --> Generates mosaic and germline scores.
                call_sv --> Based on score cut-off generates sample level files/plots for mosaic and 
                        germline mutations 
                matrix_sv --> Plot the mutation matrix
                ''')
        parser.add_argument('command', help='Please specify score/call/matrix as a command')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)
        getattr(self, args.command)()

    def cmdline(self, command):
        """Executes the command and returns the result"""
        process = Popen(args=command, stdout=PIPE, shell=True)
        return process.communicate()[0].decode()

    def comp(self, codon):
        """Returns the complement of a sequence"""
        new_cod = ""
        nuc = ["A", "G", "C", "T", "N"]
        rev_nuc = ["T", "C", "G", "A", "N"]
        for i in codon:
            new_cod = new_cod + rev_nuc[nuc.index(i.upper())]
        return new_cod

    def get_score_argument_parse(self):
        """Parses the command line arguments for get_score"""
        parser = argparse.ArgumentParser(description='get_score')
        parser.add_argument("-m", "--manifest_file",
                            help="Path to manifest file",
                            required=True,
                            type=Util.FileValidator)
        parser.add_argument("-o", "--output_dir",
                            help="Path to directory where results will be written",
                            required=True)
        parser.add_argument("-a", "--all_mutations",
                            help="Use this option to use all mutations in the vcf. By default only pass variants are "
                                 "used",
                            type=bool, nargs='?',
                            const=True, default=False)
        return parser

    def mutation_matrix_plot_argument_parse(self):
        parser = argparse.ArgumentParser(description='Plotting mutation matrix')
        parser.add_argument("-g", "--get_score_directory", help="Path to output folder of get_score", required=True,
                            type=Util.DirectoryValidator)
        parser.add_argument("-o", "--output_dir", help="Path to directory where results will be written",
                            required=True)
        parser.add_argument("-m", "--mutation", help="Mutation format should be chr_pos_alt_ref(eg. chr1_43504477_T_C)",
                            required=True, action='append')
        return parser

    def apply_score_argument_parse(self):
        """Parses the command line arguments for apply_score"""
        parser = argparse.ArgumentParser(description='apply_score')
        parser.add_argument("-g", "--get_score_directory", help="Path to output directory of the get_score option",
                            required=True, type=Util.DirectoryValidator)
        parser.add_argument("-r", "--reference", help="Path to reference file", required=True,
                            type=Util.FileValidator)
        parser.add_argument("-o", "--output_dir", help="Path to directory where results will be written",
                            required=True)
        parser.add_argument("-a", "--af_cutoff", help="Allele frequency cut-off for variants (default=0.35)",
                            default="0.35")
        parser.add_argument("-ms", "--mosaic_score_cutoff", help="Mosaic score cut-off (default=0.75)", default="0.75")
        parser.add_argument("-gs", "--germline_score_cutoff", help="Germline score cut-off (default=0.75)", default="0.75")
        parser.add_argument("-msg", "--mosaic_score_cutoff_for_germline", help="Mosaic score cut-off for germline mutations (default=0.2)", default="0.2")
        parser.add_argument("-gsm", "--germline_score_cutoff_for_mosaic", help="Germline score cut-off for mosaic mutation (default=0.5)", default="0.5")
        return parser

    def matrix(self):
        parser = self.mutation_matrix_plot_argument_parse()
        arg = parser.parse_args(sys.argv[2:])
        import seaborn as sns
        # Assigning values to variable
        ALL2_output = arg.get_score_directory
        output_dir = arg.output_dir
        mutation_list = arg.mutation
        explanation_file = os.path.join(ALL2_output, "explanation_score.txt")
        Util.ensure_dir(output_dir)

        explanation_dict = {}
        head = {}
        for i in open(explanation_file):
            line = i.strip().split("\t")
            if i.startswith("#"):
                for n, j in enumerate(line):
                    head[j] = n
                continue
            chrm = line[head["#Chrom"]]
            pos = line[head["Pos"]]
            ref = line[head["Ref"]]
            alt = line[head["Alt"]]
            mosaic_score = line[head["Mosaic_score"]]
            germline_score = line[head["Germline_score"]]
            samples = line[head["Samples_with_mutation"]].split(",")
            vaf_samples = line[head["VAF_of_samples_with_mutation"]].split(",")
            excluded_samples = line[head["Excluded_samples"]].split(",")
            mutation = "_".join([chrm, pos, ref, alt])
            mutation_related_info = {"mosaic_score": mosaic_score, "germline_score": germline_score,
                                     "sample": samples, "vaf_samples": vaf_samples, "excluded_samples":excluded_samples}
            explanation_dict[mutation] = mutation_related_info

        # loading pickle file
        print("Loading pickle file")
        mutation_matrix_file = os.path.join(ALL2_output, "mutation_matrix.pkl")
        mutation_matrix_file_fh = open(mutation_matrix_file, "rb")
        mutation_matrix_dict = pickle.load(mutation_matrix_file_fh)

        # plotting mutation
        print("Plotting mutation matrix")
        for mutation in mutation_list:
            print(mutation)
            if mutation not in mutation_matrix_dict:
                print("Mutation " + mutation + " not found")
                continue
            mutation_df = mutation_matrix_dict[mutation]
            fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [8, 1]})
            sns.heatmap(mutation_df, cmap="Blues", cbar=False, ax=ax1, linewidths=.5)
            ax1.set_ylim(len(mutation_df.index), 0)

            list_of_samples = explanation_dict[mutation]["sample"]
            list_of_samples_vaf = explanation_dict[mutation]["vaf_samples"]
            vaf_bar_list = {}
            vaf_bar_empty = {}
            for i in mutation_df.columns:
                if i in list_of_samples:
                    vaf = float(list_of_samples_vaf[list_of_samples.index(i)])
                else:
                    vaf = 0.0
                vaf_bar_list[i] = vaf
                vaf_bar_empty[i] = 1.0
            vaf_df = pd.DataFrame.from_dict(vaf_bar_list, orient='index')
            vaf_df_empty = pd.DataFrame.from_dict(vaf_bar_empty, orient='index')
            vaf_df.plot(kind='barh', xlim=(0.0, 1), legend=False, sharex=True, ax=ax2, color='g', alpha=0.5,
                        title="VAF", width=0.7)
            for n, i in enumerate(ax2.patches):
                # get_x pulls left or right; get_height pushes up or down
                ax2.text(0.3, n + 0.30, str(round(i.get_width() * 100)) + "%", fontsize=6, color='blue')

            vaf_df_empty.plot(kind='barh', xlim=(0.0, 1), legend=False,
                              sharex=True, ax=ax2, color='none', width=0.7,
                              edgecolor='blue', alpha=0.5)
            ax2.invert_yaxis()
            ax2.axis("off")
            ax2.set_title("VAF", fontsize=10)
            mosaic_score = explanation_dict[mutation]["mosaic_score"]
            germline_score = explanation_dict[mutation]["germline_score"]
            ax1.set_title(mutation+"\nMosaic_score="+str(mosaic_score)[:4]+";Germline_score="+str(germline_score)[:4])
            ## graying out unused samples
            ylabel = ax1.get_yticklabels()
            xlabel = ax1.get_xticklabels()
            excluded_samples = explanation_dict[mutation]["excluded_samples"]
            for n_x, x in enumerate(xlabel):
                if x.get_text() in excluded_samples:
                    ax1.add_patch(
                        Rectangle((n_x, 0), 1, len(ylabel), alpha=0.5, edgecolor=None, color='gray', lw=None, ls=None))
            for n_y, y in enumerate(ylabel):
                if y.get_text() in excluded_samples:
                    ax1.add_patch(
                        Rectangle((0, n_y), len(xlabel), 1, alpha=0.5, edgecolor=None, color='gray', lw=None, ls=None))
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, mutation + ".png"))
            plt.close()

    def vcf_check_for_chr(self,filename):

        if filename.endswith("vcf.gz"):
            filename_fh = gzip.open(filename)
        elif filename.endswith("vcf"):
            filename_fh = open(filename)
        for variant in filename_fh:
            try:
                variant = variant.decode("utf-8")
            except AttributeError:
                pass
            if variant.startswith("#"):
                continue
            line = variant.strip().split("\t")
            if line[0].startswith("chr"):
                return "withchr"
            else:
                return "withoutchr"

    def extract_mutation_information(self, manifest_file, output_dir, all_mutations):
        variant_dict = {}
        variant_dict_excluded ={}
        pairs_vaf_dict = {}
        head = {}
        for i in open(manifest_file):
            line = i.strip().split("\t")
            print(i)
            if i.startswith("#"):
                for n, j in enumerate(line):
                    head[j.replace("#", "")] = n
                continue
            case = line[head["Case"]]
            control = line[head["Control"]]
            try:
                case_in_vcf = line[head["Case_in_vcf"]]
                control_in_vcf = line[head["Control_in_vcf"]]
            except KeyError:
                case_in_vcf = case
                control_in_vcf = control

            try:
                bedfile = line[head["Inclusion_region"]]
            except:
                bedfile = None

            filename = line[head["Filename"]]
            pair = (case, control)
            if filename.endswith("vcf.gz"):
                filename_fh = gzip.open(filename)
            elif filename.endswith("vcf"):
                filename_fh = open(filename)
            else:
                print("Please make sure you provide a valid vcf file")
                exit()
            bed_dict={}
            if bedfile != None:
                vcf_check = self.vcf_check_for_chr(filename)
                for n,j in enumerate(open(bedfile)):
                    if j.strip() == "":
                        continue
                    bed_line = j.strip().split("\t")
                    if bed_line[0].startswith("chr") and vcf_check == "withoutchr":
                        bed_line[0]=bed_line[0].replace("chr","")
                    elif not bed_line[0].startswith("chr") and vcf_check == "withchr":
                        bed_line[0] = "chr"+bed_line[0]
                    chrm = bed_line[0]
                    start = bed_line[1]
                    end = bed_line[2]
                    if chrm in bed_dict:
                        bed_dict[chrm].append((int(start), int(end)))
                    else:
                        bed_dict[chrm] = [(int(start), int(end))]

            variant_head = {}
            # Reading variants from the vcf file
            for n,variant in enumerate(filename_fh):
                try:
                    variant = variant.decode("utf-8")
                except AttributeError:
                    pass
                line = variant.strip().split("\t")
                if variant.startswith("#"):
                    if variant.startswith("#CHROM"):
                        for n, j in enumerate(line):
                            variant_head[j] = n
                        if case_in_vcf not in variant_head or control_in_vcf not in variant_head:
                                print("Add columns to the manifest file, 'Control_in_vcf' and 'Case_in_vcf' if not already added.")
                                print("Please make sure the name of case and control match the names in the vcf file.")
                                exit()
                    continue

                filter = line[variant_head["FILTER"]]
                if filter != "PASS" and all_mutations:
                    continue
                chrm = line[variant_head["#CHROM"]]
                pos = line[variant_head["POS"]]
                ref = line[variant_head["REF"]]
                all_alt = line[variant_head["ALT"]]

                # Getting AD and DP field for case
                case_format = line[variant_head["FORMAT"]].split(":")
                try:
                    case_genotype = line[variant_head[case_in_vcf]].split(":")
                except KeyError:
                    print("Please make sure the name of the case and control in the manifest file match the"
                          " case and control specified in the vcf")
                    exit()
                try:
                    ad_index = case_format.index("AD")
                except ValueError:
                    ad_index = -1
                if ad_index > -1  and case_genotype[ad_index] != ".":
                    case_genotype_ad = case_genotype[ad_index]
                    case_genotype_depth = sum( map(int, case_genotype_ad.split(",")))
                else:
                    case_genotype_ad = "Absent"
                    case_genotype_depth = "Absent"

                if case_genotype_depth == 0:
                    continue
                # making sure to take care of multiallelic locations
                for n, alt in enumerate(all_alt.split(",")):
                    mutation = "\t".join([chrm, pos, ref, alt])
                    # storing variants
                    if mutation in variant_dict:
                        variant_dict[mutation].append(pair)
                    else:
                        variant_dict[mutation] = [pair]
                    # storing variants for excluded pairs:
                    if bedfile == None:
                        include_variant = "YES"
                    else:
                        include_variant = "NO"

                    if chrm in bed_dict:
                        for start, stop in bed_dict[chrm]:
                            if int(pos) > start and int(pos) < stop:
                                include_variant = "YES"
                                break
                    # storing VAFs
                    if case_genotype_ad != "Absent":
                        if case_genotype_depth != "Absent":
                            alt_supporting_read = case_genotype_ad.split(",")[n + 1]
                            vaf = str(float(float(alt_supporting_read) / float(case_genotype_depth)))
                        else:
                            vaf = "Absent"
                    else:
                        vaf = "Absent"

                    if include_variant == "NO":
                        if vaf == "Absent" or float(vaf) < 0.5:
                            if mutation in variant_dict_excluded:
                                variant_dict_excluded[mutation].append(pair)
                            else:
                                variant_dict_excluded[mutation] = [pair]
                    # creating a dictionary of dictionary ({pairs:{mutation:[vaf]}})
                    if pair in pairs_vaf_dict.keys():
                        if mutation in pairs_vaf_dict[pair].keys():
                            pairs_vaf_dict[pair][mutation].append(vaf)
                        else:
                            pairs_vaf_dict[pair][mutation] = [vaf]
                    else:
                        pairs_vaf_dict[pair] = {mutation: [vaf]}
        filename_fh.close()
        list_of_samples = []
        for pairs in pairs_vaf_dict:
            if pairs[0] not in list_of_samples:
                list_of_samples.append(pairs[0])
        return variant_dict, variant_dict_excluded, pairs_vaf_dict, list_of_samples

    def explanation_score(self, variant_dict, variant_dict_excluded, pairs_vaf_dict, master_list_of_samples, output_dir):

        output_file = os.path.join(output_dir, "explanation_score.txt")
        output_file_fh = open(output_file, 'w')
        output_file_fh.write("#Chrom\tPos\tRef\tAlt\tMosaic_score\tGermline_score\tTotal_samples"
                             "\tNumber_of_samples_with_mutation"
                             "\tSamples_with_mutation\tVAF_of_samples_with_mutation\tNumber_of_comparision_per_sample"
                             "\tExcluded_samples"
                             "\n")

        # number_of_cells_N is the total number of cells in the experiment
        total_number_of_cells_in_study = len(master_list_of_samples)

        # Preparing to store the data matrix for each mutation
        mutation_matrix_dict = {}
        mutation_matrix_file = os.path.join(output_dir, "mutation_matrix.pkl")
        mutation_matrix_file_fh = open(mutation_matrix_file, 'wb')

        for mutation in variant_dict:
            list_of_samples = master_list_of_samples[:]

            # pairs_all_list is a list of all pairs the mutation was called in
            pairs_all_list = variant_dict[mutation]
            # pairs_list_n is the number of all pairs the mutation was called in
            pairs_all_list_n = len(pairs_all_list)
            # pairs_excluded_list is a list of pairs the mutation was called in exclusion region
            try:
                pairs_excluded_list = variant_dict_excluded[mutation]
            except:
                pairs_excluded_list = []
            excluded_cases = []
            for case,control in pairs_excluded_list:
                if case not in excluded_cases:
                    excluded_cases.append(case)
                    if case in list_of_samples:
                        list_of_samples.remove(case)
            excluded_cases_n = len(excluded_cases)
            # pairs_excluded_list_n is the number of pairs the mutation was called in exclusion region
            pairs_excluded_list_n = len(pairs_excluded_list)

            total_number_of_cells_N = total_number_of_cells_in_study - excluded_cases_n
            pairs_list_n = pairs_all_list_n - pairs_excluded_list_n

            ## not sure waht the following line is. I think i wrote it for handling some error . Need to check later
            max_pairs_list_n = int(total_number_of_cells_N/2)*(total_number_of_cells_N-(int(total_number_of_cells_N/2)))

            # cell_fraction_f is the fraction of cells carrying the mutation

            try:
                cell_fraction_f = float(1 / 2 - float(sqrt(1 / 4 - float(pairs_list_n / (total_number_of_cells_N ** 2)))))
            except ValueError:
                print("Cell_fraction_error")
                continue
            if cell_fraction_f == 0.0:
                continue
            # is the number of cells carrying the mutation
            cells_carrying_mutation_Nv = round(cell_fraction_f * total_number_of_cells_N)
            # Creating an 'zero' data frame/matrix and updating mutation specific dataframe/matrix
            mutation_df = pd.DataFrame(np.zeros((total_number_of_cells_N, total_number_of_cells_N)),
                                       index=list_of_samples, columns=list_of_samples)
            mutation_df_for_matrix = pd.DataFrame(np.zeros((total_number_of_cells_in_study,total_number_of_cells_in_study)),
                                       index=master_list_of_samples, columns=master_list_of_samples)
            list_of_cases_with_mutation = []
            list_of_vaf_cases_with_mutation = []
            list_of_comparision_for_case = []
            list_of_excluded_cases = []
            case_dict = {}
            case_excluded_dict = {}

            for case, control in pairs_excluded_list:
                vaf = str(pairs_vaf_dict[(case, control)][mutation][0])
                if case in case_dict:
                    case_excluded_dict[case][0] = vaf
                    case_excluded_dict[case][1] = str(int(case_excluded_dict[case][1]) + 1)
                else:
                    case_excluded_dict[case] = [vaf, "1"]

            for case, control in pairs_all_list:
                vaf = str(pairs_vaf_dict[(case, control)][mutation][0])
                if case in case_dict:
                    case_dict[case][0] = vaf
                    case_dict[case][1] = str(int(case_dict[case][1]) + 1)
                else:
                    case_dict[case] = [vaf, "1"]
                mutation_df_for_matrix.loc[case, control] = 1
                if case not in case_excluded_dict and control not in case_excluded_dict:
                    mutation_df.loc[case, control] = 1

            for case in case_dict:
                list_of_cases_with_mutation.append(case)
                list_of_vaf_cases_with_mutation.append(case_dict[case][0])
                if case not in case_excluded_dict:
                    list_of_comparision_for_case.append(case_dict[case][1])
            for case in case_excluded_dict:
                list_of_excluded_cases.append(case)
            if list_of_excluded_cases == []:
                excluded_cases = "-"
            else:
                excluded_cases = ",".join(list_of_excluded_cases)

            mutation_matrix_dict["_".join(mutation.split("\t"))] = mutation_df_for_matrix

            # calculating explanation score

            ordered_col_sum = mutation_df.sum(axis=1).sort_values(ascending=False)


            ordered_row_sum = mutation_df.sum(axis=0).sort_values(ascending=False)

            explained_call_n_mosaic = ordered_col_sum[:cells_carrying_mutation_Nv].sum()
            explained_call_n_germ = ordered_row_sum[:cells_carrying_mutation_Nv].sum()
            explanation_score_mosaic = explained_call_n_mosaic / pairs_list_n
            explanation_score_germ = explained_call_n_germ / pairs_list_n
            cells_with_mutation = ",".join(list_of_cases_with_mutation)
            vaf_of_samples_with_mutation = ",".join(list_of_vaf_cases_with_mutation)
            comparision_for_case = ",".join(list_of_comparision_for_case)
            if cells_with_mutation == "":
                cells_with_mutation = "-"
            output_line = "\t".join(["\t".join(mutation.split("\t")),
                                     str(explanation_score_mosaic),
                                     str(explanation_score_germ),
                                     str(total_number_of_cells_in_study),
                                     str(len(list_of_cases_with_mutation)),
                                     cells_with_mutation,
                                     vaf_of_samples_with_mutation,
                                     comparision_for_case,
                                     excluded_cases])
            output_file_fh.write(output_line + "\n")

        pickle.dump(mutation_matrix_dict, mutation_matrix_file_fh)
        mutation_matrix_file_fh.close()
        output_file_fh.close()

    def plot_score(self, output_dir):
        explanation_score = os.path.join(output_dir, "explanation_score.txt")

        # Plotting germline versus mosaic score scatter plot
        df_es = pd.read_csv(explanation_score, sep="\t")
        x = df_es["Mosaic_score"]
        y = df_es["Germline_score"]
        size_dict = {}
        for n, i in enumerate(x):
            score_pair = str(i) + "_" + str(y[n])
            if score_pair in size_dict:
                size_dict[score_pair] += 1
            else:
                size_dict[score_pair] = 1
        size = []
        for n, i in enumerate(x):
            score_pair = str(i) + "_" + str(y[n])
            size.append(float(size_dict[score_pair]))
        df_es.plot.scatter("Mosaic_score", "Germline_score", c=size, s=[float(s / 10.0) for s in size], cmap='tab10', colorbar = True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "Explanation_score_scatter.png"))
        plt.close()

    def plot_curve(self, x, mosaic_score_cut, germ_score_cut):

        y_curve = np.multiply(np.sqrt(np.subtract(1,
                                       np.divide(np.power(x, 2),
                                        np.power(mosaic_score_cut, 2)))),
                              germ_score_cut)


        return y_curve

    def plot_score_annotate(self, explanation_score_file, output_dir, mosaic_score, germline_score,mosaic_cutoff_for_germline_mutations,
                                 germline_cutoff_for_mosaic_mutations):
        explanation_score = explanation_score_file
        mosaic_score_cut = mosaic_score
        germ_score_cut = germline_score
        germ_score_mosaic_cut = mosaic_cutoff_for_germline_mutations    # Vivek needs to change this to 90 percentile
        mosaic_score_germ_cut = germline_cutoff_for_mosaic_mutations    # Vivek needs to change this to 90 percentile
        # Plotting germline versus mosaic score scatter plot
        df_es = pd.read_table(explanation_score)
        x = df_es["Mosaic_score"]
        y = df_es["Germline_score"]
        size_dict = {}
        points_to_plot = 1000
        for n, i in enumerate(x):
            score_pair = str(i) + "_" + str(y[n])
            if score_pair in size_dict:
                size_dict[score_pair] += 1
            else:
                size_dict[score_pair] = 1
        size = []
        for n, i in enumerate(x):
            score_pair = str(i) + "_" + str(y[n])
            size.append(float(size_dict[score_pair]))
        ax = df_es.plot.scatter("Mosaic_score", "Germline_score", c=size, s=[float(s / 10.0) for s in size], cmap = 'tab10', colorbar = True)
        x_curve = np.linspace(0.0, mosaic_score_cut, points_to_plot)
        y_curve = self.plot_curve(x_curve, mosaic_score_cut, germ_score_cut)
        line_color = 'black'
        # Noise
        ax.fill_between(x_curve, y_curve, color='r', alpha=0.1)

        ## Germline
        x_germ = np.linspace(0.0, germ_score_mosaic_cut, points_to_plot)
        y_germ = self.plot_curve(x_germ, mosaic_score_cut, germ_score_cut)
        y2_germ = np.linspace(1, 1, points_to_plot)
        ax.fill_between(x_germ, y_germ, y2_germ, color='b', alpha=0.1)

        ## high freq Mosaic
        x_end_point = self.plot_curve([mosaic_score_germ_cut], germ_score_cut, mosaic_score_cut)[0]  # 0.5
        x_mosaic_high = np.linspace(germ_score_mosaic_cut, 1, points_to_plot)
        y_mosaic_high = []
        for i in x_mosaic_high:
            if i > x_end_point:
                y_mosaic_high.append(mosaic_score_germ_cut)
            else:
                y_mosaic_high.append(self.plot_curve([i], mosaic_score_cut, germ_score_cut)[0])
        y2_mosaic_high = np.linspace(1, 1, points_to_plot)
        ax.fill_between(x_mosaic_high, y_mosaic_high, y2_mosaic_high, color='black', alpha=0.1)

        # Mosaic
        x_start_point = self.plot_curve([mosaic_score_germ_cut], germ_score_cut, mosaic_score_cut)[0]
        x_end_point = self.plot_curve([0], germ_score_cut, mosaic_score_cut)[0]
        x_mosaic = np.linspace(x_start_point, 1, points_to_plot)
        y_mosaic = self.plot_curve(x_mosaic, mosaic_score_cut, germ_score_cut)
        y_mosaic = []
        for i in x_mosaic:
            if i < x_end_point:
                y_mosaic.append(self.plot_curve([i], mosaic_score_cut, germ_score_cut)[0])
            elif i >= x_end_point:
                y_mosaic.append(0.0)
        y2_mosaic = np.linspace(mosaic_score_germ_cut, mosaic_score_germ_cut, points_to_plot)
        ax.fill_between(x_mosaic, y_mosaic, y2_mosaic, color='green', alpha=0.1)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "Explanation_score_scatter_annotated.png"))
        plt.close()

    def score(self):
        # Extracting passed arguments
        parser = self.get_score_argument_parse()
        arg = parser.parse_args(sys.argv[2:])

        # Assigning values to variable
        manifest_file = arg.manifest_file
        output_dir = arg.output_dir
        all_mutations = arg.all_mutations
        Util.ensure_dir(output_dir)

        # Extracting variant information from the manifest file.
        print("Extracting variant information")
        variant_dict, variant_dict_excluded, pairs_vaf_dict, list_of_samples = self.extract_mutation_information(manifest_file, output_dir,all_mutations)
        # variant_dict={mutation:[(case,control)]}
        # pairs_vaf_dict={pairs:{mutation:[vaf]}}
        # list_of_samples = list of all samples in the analysis

        # Generating explanation score
        print("Generating explanation scores")
        self.explanation_score(variant_dict, variant_dict_excluded, pairs_vaf_dict, list_of_samples, output_dir)

        # Plotting
        print("Plotting")
        # plotting(output_dir,list_of_samples,af_cutoff)
        self.plot_score(output_dir)

    def plot_bar(self, vaf_dict, af_cutoff, output_dir):
        """ plotting bar plot """
        output_file = os.path.join(output_dir, "mutation_type_count.png")
        mutation_count = {"Mosaic": 0, "Germline": 0, "Noise": 0}
        mutation_count_per_sample = {}
        variant_list = []
        for sample in vaf_dict:
            mutation_count_per_sample[sample] = {"Mosaic": 0, "Germline": 0, "Noise": 0}
            for variant_type in ["Mosaic", "Germline", "Noise"]:
                for pos, vaf, germline_score, mosaic_score in vaf_dict[sample][variant_type]:
                    if pos not in variant_list:
                        mutation_count[variant_type] += 1
                        variant_list.append(pos)
                    if vaf >= af_cutoff:
                        mutation_count_per_sample[sample][variant_type] += 1
        df = pd.DataFrame.from_dict(mutation_count, orient='index')
        title = "Number of unique mutations across all samples"
        ax = df.plot(kind='bar', legend=False, title=title, alpha=0.60, grid=False)

        for i in ax.patches:
            # get_x pulls left or right; get_height pushes up or down
            ax.text(i.get_x() + 0.15, i.get_height() / 2,
                    str(i.get_height()), fontsize=10, color='black')

        plt.ylabel("Number of mutations")
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()

        # Per sample variant counts
        df = pd.DataFrame.from_dict(mutation_count_per_sample, orient='index')
        output_file_per_sample = os.path.join(output_dir, "per_sample_mutation_count.png")
        title = "Mutations with allele frequency cutoff above " + str(af_cutoff)
        ax = df.plot(kind='bar', legend=True, grid=False)
        plt.title(title)
        plt.ylabel("Number of mutations")
        plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.5))
        plt.tight_layout()
        plt.savefig(output_file_per_sample)
        plt.close()

    def plot_vaf(self, vaf_dict, af_cutoff, output_dir):
        """ plotting vaf plots """
        for sample in vaf_dict:
            mutation_list_snp = {"Mosaic": [], "Germline": [], "Noise": []}
            mutation_list_indel = {"Mosaic": [], "Germline": [], "Noise": []}
            for variant_type in ["Mosaic", "Germline", "Noise"]:
                output_file = os.path.join(output_dir, sample + "." + variant_type + ".vaf_plot.png")
                for pos, vaf, germline_score, mosaic_score in vaf_dict[sample][variant_type]:
                    chrm, position, ref, alt = pos.split("_")
                    if len(ref) > 1 or len(alt) > 1:
                        mutation_list_indel[variant_type].append(vaf)
                    else:
                        mutation_list_snp[variant_type].append(vaf)

                title = variant_type + " VAF distribution for " + sample
                df_snp = pd.DataFrame.from_dict(mutation_list_snp[variant_type])
                df_indel = pd.DataFrame.from_dict(mutation_list_indel[variant_type])
                if df_snp.empty == False and df_indel.empty == False:
                    ax1 = df_snp.plot(kind='hist', bins=33, range=(0, 1),
                                alpha = 1, legend=None,
                                weights=np.ones_like(df_snp[df_snp.columns[0]]) * 1. / len(df_snp),
                                title = title)
                    ax1 = df_indel.plot(kind='hist', bins=33, range=(0, 1),
                                  alpha=0.5, legend=None,
                                  weights=np.ones_like(df_indel[df_indel.columns[0]]) * 1. / len(df_indel),
                                  ax = ax1)

                    plt.legend(["SNV", "INDEL"])
                elif df_indel.empty:
                    df_snp.plot(kind='hist', bins=33, range=(0, 1),
                                alpha=0.5, legend=None,
                                weights=np.ones_like(df_snp[df_snp.columns[0]]) * 1. / len(df_snp),
                                title=title)
                    plt.legend(["SNV"])
                elif df_snp.empty:
                    df_indel.plot(kind='hist', bins=33, range=(0, 1),
                                  alpha=0.5, legend=None,
                                  weights=np.ones_like(df_indel[df_indel.columns[0]]) * 1. / len(df_indel),
                                  title=title)
                    plt.legend(["INDEL"])
                plt.ylim(0.0, 0.3)
                plt.ylabel("Percent of mutations")
                plt.axvline(0.5, color='r', linestyle='dashed', linewidth=1)
                plt.axvline(af_cutoff, color='g', linestyle='dashed', linewidth=1)
                plt.tight_layout()
                plt.savefig(output_file)
                plt.close()

    def create_spectra_file(self, vaf_dict, sample, variant_type, af_cutoff, output_dir, reference):
        """generates the underlying file for spectrum plot"""
        sample_signature_count = {}
        total_variant = 0
        for pos, vaf, germline_score, mosaic_score in vaf_dict[sample][variant_type]:
            line = pos.strip().split("_")
            chrm = line[0]
            pos = line[1]
            ref = line[2].upper()
            alt = line[3].upper()
            if len(ref) > 1 or len(alt) > 1 or vaf < af_cutoff:
                continue
            ref_extract_pos = chrm + ":" + str(int(pos) - 1) + "-" + str(int(pos) + 1)
            cmd = " ".join(["samtools faidx", reference,
                            ref_extract_pos,
                            "|grep -v \"^>\""])

            tri_base = self.cmdline(cmd).strip().upper()
            if ref == "G" or ref == "A":
                ref = self.comp(ref)
                alt = self.comp(alt)
                tri_base = self.comp(tri_base)
            tri_base_sig = tri_base + ":" + alt
            if tri_base_sig in sample_signature_count:
                sample_signature_count[tri_base_sig] += 1
            else:
                sample_signature_count[tri_base_sig] = 1
            total_variant += 1
        mutation_spectra_file = os.path.join(output_dir, sample + "." + variant_type + ".mutation_spectra.txt")
        mutation_spectra_file_fh = open(mutation_spectra_file, 'w')
        # print total_variant
        script_path = os.path.dirname(sys.argv[0])
        list_of_mut = os.path.join(script_path, "list_of_tri_bases.txt")
        for i in open(list_of_mut):
            if i.startswith("#"):
                continue
            mut_sig = i.strip()
            if mut_sig in sample_signature_count:
                spectra_line = mut_sig + "\t" + str(
                    float(float(sample_signature_count[mut_sig]) / float(total_variant)) * 100) + "\t" + str(
                    sample_signature_count[mut_sig]) + "\n"
            else:
                spectra_line = mut_sig + "\t0.0\t0\n"
            mutation_spectra_file_fh.write(spectra_line)
        return mutation_spectra_file

    def six_mutation_spectrum_plot(self, mutation_spectra_file, output_dir, sample_name, variant_type):
        six_mutation_file = os.path.join(output_dir, sample_name + "." + variant_type + ".six_mutation_spectrum.txt")
        six_mutation_plot = os.path.join(output_dir, sample_name + "." + variant_type + ".six_mutation_spectrum.png")
        six_mutations = {}
        six_mutation_count = {}
        cpg = 0
        cpg_per = 0.0
        import matplotlib.patches as mpatches
        for i in open(mutation_spectra_file):
            line = i.strip().split("\t")
            ref = line[0][1]
            ref_plus_one = line[0][2]
            alt = line[0][4]
            mutation = ref + ">" + alt
            per = float(line[1])
            count = int(line[2])
            if mutation in six_mutations:
                if mutation == "C>T" and ref_plus_one == "G":
                    cpg += count
                    cpg_per += per
                else:
                    six_mutations[mutation] += per
                    six_mutation_count[mutation] += count
            else:
                six_mutations[mutation] = per
                six_mutation_count[mutation] = count
        six_mutation_fh = open(six_mutation_file, 'w')
        labels = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
        six_mutation_fh.write("Mutation\tMutation_percent\tMutation_count\tCpG_percent\tCpG_count\n")
        for i in labels:
            if i == "C>T":
                out = i + "\t" + str(six_mutations[i]) + "\t" + str(six_mutation_count[i]) + "\t" + str(
                    cpg_per) + "\t" + str(cpg) + "\n"
            else:
                out = i + "\t" + str(six_mutations[i]) + "\t" + str(six_mutation_count[i]) + "\t0.0\t0\n"
            six_mutation_fh.write(out)

        mycolor = ["deepskyblue", "black", "red", "silver", "lightgreen", "lightpink"]

        six_mutation_fh.close()
        ylim_max = 90
        df = pd.read_csv(six_mutation_file, sep="\t", header=0)
        df["All_mutation_per"] = df["Mutation_percent"] + df["CpG_percent"]
        ax = df[["Mutation", "All_mutation_per"]].plot(kind='bar', grid=0, legend=False, color="red", alpha=0.8,
                                                       hatch="////", figsize=(15, 3), ylim=(0.0, ylim_max))
        df[["Mutation", "Mutation_percent"]].plot(kind='bar', grid=0, legend=False, color=[mycolor], figsize=(15, 3),
                                                  ylim=(0.0, ylim_max), ax=ax)
        a = ax.set_xticklabels(labels, size=10)
        plt.ylabel("Percent of mutation")
        plt.title(variant_type + " six Mutation Signature for " + sample_name, y=1.10)
        cpg_patch = mpatches.Patch(facecolor='red', label='CpG', hatch="////")
        plt.legend(handles=[cpg_patch], bbox_to_anchor=(0.53, 0.60))
        plt.tight_layout()
        plt.savefig(six_mutation_plot)
        plt.close()

    def mutation_spectrum_plot(self, mutation_spectra_file, output_dir, sample_name, variant_type):
        """Creates the spectrum plot"""
        global color_sig
        ylim_max = 15.0
        df = pd.read_table(mutation_spectra_file, header=None)
        df = df[df.columns[0:2]]
        label = []
        top_label = []
        mycolor = []
        for i in df.iloc[:, 0]:
            mut_sig = i.strip().split("\t")[0]
            ref = mut_sig[1].upper()
            alt = mut_sig[4].upper()
            label.append(mut_sig.split(":")[0])
            top_label.append(ref + ">" + alt)
            if ref == "C" and alt == "A":
                color_sig = "deepskyblue"
            if ref == "C" and alt == "G":
                color_sig = "black"
            if ref == "C" and alt == "T":
                color_sig = "red"
            if ref == "T" and alt == "A":
                color_sig = "silver"
            if ref == "T" and alt == "C":
                color_sig = "lightgreen"
            if ref == "T" and alt == "G":
                color_sig = "lightpink"

            mycolor.append(color_sig)

        ax = df.plot(kind='bar', grid=0, legend=False, color=[mycolor], ylim=(0.0, ylim_max), figsize=(15, 3))
        plt.plot([0, 13.5], [ylim_max, ylim_max], 'k-', lw=35, color='deepskyblue')
        plt.title(variant_type + " mutational Signature for " + sample_name, y=1.10)
        plt.ylabel("Percent of mutation")
        plt.text(6, ylim_max + 1, "C>A", fontsize=12)
        plt.plot([17, 29.5], [ylim_max, ylim_max], 'k-', lw=35, color='black')
        plt.text(22, ylim_max + 1, "C>G", fontsize=12)
        plt.plot([33, 46], [ylim_max, ylim_max], 'k-', lw=35, color='red')
        plt.text(38, ylim_max + 1, "C>T", fontsize=12)

        plt.plot([49.5, 62], [ylim_max, ylim_max], 'k-', lw=35, color='silver')
        plt.text(54, ylim_max + 1, "T>A", fontsize=12)
        plt.plot([65.5, 78], [ylim_max, ylim_max], 'k-', lw=35, color='lightgreen')
        plt.text(70, ylim_max + 1, "T>C", fontsize=12)
        plt.plot([81.5, 95], [ylim_max, ylim_max], 'k-', lw=35, color='lightpink')
        plt.text(85, ylim_max + 1, "T>G", fontsize=12)
        output_plot = os.path.join(output_dir, sample_name + "." + variant_type + ".mutation_spectrum.png")
        a = ax.set_xticklabels(label, size=10)
        plt.tight_layout()
        plt.savefig(output_plot)
        plt.close()

    def plot_mutation_spectrum(self, vaf_dict, af_cutoff, output_dir, reference):
        """ Plots mutation spectrum """
        for sample in vaf_dict:
            mutation_list = {"Mosaic": [], "Germline": [], "Noise": []}
            for variant_type in ["Mosaic", "Germline", "Noise"]:
                output_file = os.path.join(output_dir, sample + "." + variant_type + ".vaf_plot.png")
                mutation_spectra_file = self.create_spectra_file(vaf_dict, sample, variant_type, af_cutoff,
                                                                 output_dir, reference)
                self.mutation_spectrum_plot(mutation_spectra_file, output_dir, sample, variant_type)
                self.six_mutation_spectrum_plot(mutation_spectra_file, output_dir, sample, variant_type)

    def per_sample_mutation(self, vaf_dict, per_sample_mutation_dir, af_cutoff):
        """Creating per sample mutation file"""
        for sample in vaf_dict:
            file_out = os.path.join(per_sample_mutation_dir, sample+".tsv")
            file_out_fh = open(file_out,'w')
            file_out_fh.write("#Chr\tPos\tRef\tAlt\tVAF\tMosaic_score\tGermline_score\tVariant_type\n")
            for variant_type in vaf_dict[sample]:
                for pos, vaf, germline_score, mosaic_score in vaf_dict[sample][variant_type]:
                    if vaf < af_cutoff:
                        continue
                    file_out_fh.write("\t".join(["\t".join(pos.split("_")), str(vaf), str(mosaic_score), str(germline_score), variant_type])+"\n")
            file_out_fh.close()

    def call(self):
        # Extracting passed arguments
        parser = self.apply_score_argument_parse()
        arg = parser.parse_args(sys.argv[2:])

        # Assigning values to variable
        get_score_dir = arg.get_score_directory
        output_dir = arg.output_dir
        reference = arg.reference
        af_cutoff = float(arg.af_cutoff)
        mosaic_score_cutoff = float(arg.mosaic_score_cutoff)
        germline_score_cutoff = float(arg.germline_score_cutoff)


        Util.ensure_dir(output_dir)
        explanation_score_file = os.path.join(get_score_dir, "explanation_score.txt")
        vaf_dict = {}
        head = {}

        df_manifest = pd.read_table(explanation_score_file)

        mosaic_cutoff_for_germline_mutations = float(arg.mosaic_score_cutoff_for_germline)
        germline_cutoff_for_mosaic_mutations = float(arg.germline_score_cutoff_for_mosaic)

        for i in open(explanation_score_file):
            line = i.strip().split("\t")
            if i.startswith("#"):
                for n, j in enumerate(line):
                    head[j] = n
                continue
            pos = "_".join(line[head["#Chrom"]:head["Mosaic_score"]])
            list_of_samples = line[head["Samples_with_mutation"]].split(",")
            list_of_vafs = line[head["VAF_of_samples_with_mutation"]].split(",")
            variant_type = ""
            germline_score = float(line[head["Germline_score"]])
            mosaic_score = float(line[head["Mosaic_score"]])

            # this block annotated the mutation as mosaic, germline or noise
            if mosaic_score > mosaic_score_cutoff:
                y_curve = y_curve = self.plot_curve([mosaic_score_cutoff], mosaic_score_cutoff, germline_score_cutoff)[0]
            else:
                y_curve = self.plot_curve([mosaic_score], mosaic_score_cutoff, germline_score_cutoff)[0]
            if germline_score > germline_score_cutoff:
                x_curve = self.plot_curve([germline_score_cutoff], germline_score_cutoff, mosaic_score_cutoff)[0]
            else:
                x_curve = self.plot_curve([germline_score], germline_score_cutoff, mosaic_score_cutoff)[0]
            if mosaic_score < x_curve and germline_score < y_curve:
                variant_type = "Noise"
            elif germline_score > y_curve and mosaic_score > x_curve and mosaic_score >= mosaic_cutoff_for_germline_mutations and germline_score > germline_cutoff_for_mosaic_mutations:
                variant_type = "Mosaic_high_freq"
            elif germline_score > y_curve and mosaic_score > x_curve and mosaic_score >= mosaic_cutoff_for_germline_mutations and germline_score <= germline_cutoff_for_mosaic_mutations:
                variant_type = "Mosaic"
            elif germline_score > y_curve and x_curve < mosaic_score < mosaic_cutoff_for_germline_mutations:
                variant_type = "Germline"
            else:
                print("variant type not detected, check for bugs in variant type definition")

            for index, sample in enumerate(list_of_samples):
                vaf = float(list_of_vafs[index])
                if sample in vaf_dict:
                    if variant_type in vaf_dict[sample]:
                        vaf_dict[sample][variant_type].append((pos, vaf, germline_score, mosaic_score))
                    else:
                        vaf_dict[sample][variant_type] = [(pos, vaf, germline_score, mosaic_score)]
                else:
                    vaf_dict[sample] = {}

        # structure of var_dict={sample:{variant_type:[(pos,vaf)]}}  ; pos = chr pos ref alt

        # per sample mutation file
        per_sample_mutation_dir = os.path.join(output_dir, "per_sample_mutation")
        Util.ensure_dir(per_sample_mutation_dir)
        self.per_sample_mutation(vaf_dict,per_sample_mutation_dir,af_cutoff)



        # plotting
        self.plot_score_annotate(explanation_score_file, output_dir, mosaic_score_cutoff, germline_score_cutoff, mosaic_cutoff_for_germline_mutations,
                                 germline_cutoff_for_mosaic_mutations)

        mutation_count_output_dir = os.path.join(output_dir, "mutation_counts")
        Util.ensure_dir(mutation_count_output_dir)
        self.plot_bar(vaf_dict, af_cutoff, mutation_count_output_dir)

        vaf_output_dir = os.path.join(output_dir, "vaf_plots")
        Util.ensure_dir(vaf_output_dir)
        self.plot_vaf(vaf_dict, af_cutoff, vaf_output_dir)

        mutation_spectrum_output_dir = os.path.join(output_dir, "mutation_spectrum")
        Util.ensure_dir(mutation_spectrum_output_dir)
        self.plot_mutation_spectrum(vaf_dict, af_cutoff, mutation_spectrum_output_dir, reference)


    # SV module starts here

    def reciprocal_overlap(self, SV_dict, chr_start_end_svtype, overlap_percent = 0.5):
        """ 50 % reciprocal overlap"""
        chr = chr_start_end_svtype.split("\t")[0]
        start = int(chr_start_end_svtype.split("\t")[1])
        end = int(chr_start_end_svtype.split("\t")[2])
        svtype = chr_start_end_svtype.split("\t")[3]
        sv_len = end - start
        sv_list = range(start,end+1)

        for sv in SV_dict:
            for mutation in SV_dict[sv]:
                dict_chr_start_end_svtype = mutation.split("\t")
                dict_chr = dict_chr_start_end_svtype[0]
                dict_start = int(dict_chr_start_end_svtype[1])
                dict_end = int(dict_chr_start_end_svtype[2])
                dict_svtype = dict_chr_start_end_svtype[3]
                dict_sv_len = dict_end - dict_start
                if chr != dict_chr or svtype != dict_svtype:
                    break

                if start >= dict_start and start < dict_end:
                    overlap_start = start
                elif dict_start >= start and dict_start < end:
                    overlap_start = dict_start
                else:
                    continue
                if end < dict_end:
                    overlap_end = end
                else:
                    overlap_end = dict_end

                overlap = overlap_end - overlap_start

                if overlap >= int(sv_len * overlap_percent) and overlap >= int(dict_sv_len * overlap_percent):
                    return sv
        return False

    def get_score_argument_parse_sv(self):
        """Parses the command line arguments for get_score"""
        parser = argparse.ArgumentParser(description='get_score')
        parser.add_argument("-m", "--manifest_file",
                            help="Path to manifest file",
                            required=True,
                            type=Util.FileValidator)
        parser.add_argument("-o", "--output_dir",
                            help="Path to directory where results will be written",
                            required=True)
        parser.add_argument("-a", "--all_mutations",
                            help="Use this option to use all mutations in the vcf. By default only pass variants are "
                                 "used",
                            type=bool, nargs='?',
                            const=True, default=False)
        return parser

    def mutation_matrix_plot_argument_parse_sv(self):
        parser = argparse.ArgumentParser(description='Plotting mutation matrix')
        parser.add_argument("-g", "--get_score_directory", help="Path to output folder of get_score", required=True,
                            type=Util.DirectoryValidator)
        parser.add_argument("-o", "--output_dir", help="Path to directory where results will be written",
                            required=True)
        parser.add_argument("-m", "--mutation", help="Provide the SV number(eg. 27)",
                            required=True, action='append')
        return parser

    def apply_score_argument_parse_sv(self):
        """Parses the command line arguments for apply_score"""
        parser = argparse.ArgumentParser(description='apply_score')
        parser.add_argument("-g", "--get_score_directory", help="Path to output directory of the get_score option",
                            required=True, type=Util.DirectoryValidator)
        parser.add_argument("-o", "--output_dir", help="Path to directory where results will be written",
                            required=True)
        parser.add_argument("-ms", "--mosaic_score_cutoff", help="Mosaic score cut-off (default=0.75)", default="0.75")
        parser.add_argument("-gs", "--germline_score_cutoff", help="Germline score cut-off (default=0.75)", default="0.75")
        parser.add_argument("-msg", "--mosaic_score_cutoff_for_germline", help="Mosaic score cut-off for germline mutations (default=0.2)", default="0.2")
        parser.add_argument("-gsm", "--germline_score_cutoff_for_mosaic", help="Germline score cut-off for mosaic mutation (default=0.5)", default="0.5")
        return parser

    def matrix_sv(self):
        parser = self.mutation_matrix_plot_argument_parse_sv()
        arg = parser.parse_args(sys.argv[2:])
        import seaborn as sns
        # Assigning values to variable
        ALL2_output = arg.get_score_directory
        output_dir = arg.output_dir
        mutation_list = arg.mutation
        explanation_file = os.path.join(ALL2_output, "explanation_score.txt")
        Util.ensure_dir(output_dir)

        explanation_dict = {}
        head = {}
        for i in open(explanation_file):
            line = i.strip().split("\t")
            if i.startswith("#"):
                for n, j in enumerate(line):
                    head[j] = n
                continue
            mosaic_score = line[head["Mosaic_score"]]
            germline_score = line[head["Germline_score"]]
            samples = line[head["Samples_with_mutation"]].split(",")
            mutation = line[head["#SV"]]
            mutation_related_info = {"mosaic_score": mosaic_score, "germline_score": germline_score,
                                     "sample": samples}
            explanation_dict[mutation] = mutation_related_info
        # loading pickle file
        print("Loading pickle file")
        mutation_matrix_file = os.path.join(ALL2_output, "mutation_matrix.pkl")
        mutation_matrix_file_fh = open(mutation_matrix_file, "rb")
        mutation_matrix_dict = pickle.load(mutation_matrix_file_fh)

        # plotting mutation
        print("Plotting mutation matrix")
        for mutation in mutation_list:
            print(mutation)
            if mutation not in mutation_matrix_dict:
                print("Mutation " + mutation + " not found")
                continue
            mutation_df = mutation_matrix_dict[mutation]
            ax1 = sns.heatmap(mutation_df, cmap="Blues", cbar=False, linewidths=.5)
            ax1.set_ylim(len(mutation_df.index), 0)

            ax1.set_title(mutation)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, mutation + ".png"))
            plt.close()

    def extract_mutation_information_sv(self, manifest_file, output_dir, all_mutations):
        variant_dict = {}
        SV_mutations_dict = {}
        SV_dict = {}
        pairs_vaf_dict = {}
        head = {}
        sv_count = 0
        for i in open(manifest_file):
            line = i.strip().split("\t")
            if i.startswith("#"):
                for n, j in enumerate(line):
                    head[j.replace("#", "")] = n
                continue
            case = line[head["Case"]]
            control = line[head["Control"]]
            try:
                case_in_vcf = line[head["Case_in_vcf"]]
                control_in_vcf = line[head["Control_in_vcf"]]
            except KeyError:
                case_in_vcf = case
                control_in_vcf = control
            filename = line[head["Filename"]]
            pair = (case, control)
            if filename.endswith("vcf.gz"):
                filename_fh = gzip.open(filename)
            elif filename.endswith("vcf"):
                filename_fh = open(filename)
            else:
                print("Please make sure you provide a valid vcf file")
                exit()
            variant_head = {}
            # Reading variants from the vcf file
            for variant in filename_fh:
                try:
                    variant = variant.decode("utf-8")
                except AttributeError:
                    pass
                line = variant.strip().split("\t")
                if variant.startswith("#"):
                    if variant.startswith("#CHROM"):
                        for n, j in enumerate(line):
                            variant_head[j] = n
                        if case_in_vcf not in variant_head or control_in_vcf not in variant_head:
                                print("Add columns to the manifest file, 'Control_in_vcf' and 'Case_in_vcf' if not already added.")
                                print("Please make sure the name of case and control match the names in the vcf file.")
                                exit()
                    continue
                filter = line[variant_head["FILTER"]]
                if filter != "PASS" and all_mutations:
                    continue
                chrm = line[variant_head["#CHROM"]]
                start_pos = line[variant_head["POS"]]
                for j in line[variant_head["INFO"]].split(";"):
                    if j.startswith("END"):
                        end_pos = j.split("=")[1]
                    if j.startswith("SVTYPE"):
                        sv_type = j.split("=")[1]
                if sv_type == "BND":
                    continue
                chr_start_end_svtype = chrm + "\t" + start_pos + "\t" + end_pos + "\t" + sv_type
                ref = line[variant_head["REF"]]
                alt = line[variant_head["ALT"]]
                mutation = "\t".join([chrm, start_pos, ref, alt, sv_type])
                sv = self.reciprocal_overlap(SV_dict, chr_start_end_svtype)
                if sv == False:
                    sv_count += 1
                    sv = str(sv_count)
                    SV_dict[sv] = [chr_start_end_svtype]
                    SV_mutations_dict[sv] = {pair:[mutation]}
                else:
                    SV_dict[sv].append(chr_start_end_svtype)
                    if pair in SV_mutations_dict[sv]:
                        if mutation not in SV_mutations_dict[sv][pair]:
                            SV_mutations_dict[sv][pair].append(mutation)
                    else:
                        SV_mutations_dict[sv][pair] = [mutation]

                # Getting AD and DP field for case
                case_format = line[variant_head["FORMAT"]].split(":")
                try:
                    case_genotype = line[variant_head[case_in_vcf]].split(":")
                except KeyError:
                    print("Please make sure the name of the case and control in the manifest file match the"
                          " case and control specified in the vcf")
                    exit()
                if sv not in variant_dict:
                    variant_dict[sv]=[pair]
                else:
                    variant_dict[sv].append(pair)
                vaf = 0
                if pair in pairs_vaf_dict.keys():
                    if sv in pairs_vaf_dict[pair].keys():
                        pairs_vaf_dict[pair][sv].append(vaf)
                    else:
                        pairs_vaf_dict[pair][sv] = [vaf]
                else:
                    pairs_vaf_dict[pair] = {sv: [vaf]}
        filename_fh.close()
        list_of_samples = []
        for pairs in pairs_vaf_dict:
            if pairs[0] not in list_of_samples:
                list_of_samples.append(pairs[0])
        return variant_dict, SV_mutations_dict, list_of_samples

    def explanation_score_sv(self, variant_dict, SV_mutation_dict, list_of_samples, output_dir):

        output_file = os.path.join(output_dir, "explanation_score.txt")
        output_sv_mapping_file = os.path.join(output_dir, "SV_mapping.txt")
        output_file_fh = open(output_file, 'w')
        output_file_fh.write("#SV\tMosaic_score\tGermline_score\tNumber_of_samples_with_mutation"
                             "\tSamples_with_mutation\tNumber_of_comparision_per_sample\tSV_type"
                             "\n")

        output_sv_mapping_file_fh = open(output_sv_mapping_file, 'w')
        output_sv_mapping_file_fh.write("#SV\tAssociated_mutations\n")

        # number_of_cells_N is the total number of cells in the experiment
        total_number_of_cells_N = len(list_of_samples)

        # Preparing to store the data matrix for each mutation
        mutation_matrix_dict = {}
        mutation_matrix_file = os.path.join(output_dir, "mutation_matrix.pkl")
        mutation_matrix_file_fh = open(mutation_matrix_file, 'wb')

        for sv in variant_dict:
            # pairs_list is a list of pairs the mutation was called in
            pairs_list = variant_dict[sv]
            # pairs_list_n is the number of pairs the mutation was called in
            pairs_list_n = len(pairs_list)
            max_pairs_list_n = int(total_number_of_cells_N/2)*(total_number_of_cells_N-(int(total_number_of_cells_N/2)))
            # cell_fraction_f is the fraction of cells carrying the mutation
            try:
                cell_fraction_f = 1 / 2 - sqrt(1 / 4 - pairs_list_n / total_number_of_cells_N ** 2)
            except ValueError:
                continue
            # is the number of cells carrying the mutation
            cells_carrying_mutation_Nv = round(cell_fraction_f * total_number_of_cells_N)
            # Creating an 'zero' data frame/matrix and updating mutation specific dataframe/matrix
            mutation_df = pd.DataFrame(np.zeros((total_number_of_cells_N, total_number_of_cells_N)),
                                       index=list_of_samples, columns=list_of_samples)
            list_of_cases_with_mutation = []
            list_of_comparision_for_case = []
            case_dict = {}
            list_of_mutation = []
            sv_type = ""
            for case, control in pairs_list:
                for mutation in SV_mutation_dict[sv][(case, control)]:
                        mapping_value = case+":"+":".join(mutation.split("\t"))
                        if mapping_value not in list_of_mutation:
                            list_of_mutation.append(mapping_value)
                        sv_type = mutation.split("\t")[-1]
                if case in case_dict:
                    case_dict[case][0] = str(int(case_dict[case][0]) + 1)
                else:
                    case_dict[case] = ["1"]
                mutation_df.loc[case, control] = 1

            for case in case_dict:
                list_of_cases_with_mutation.append(case)
                list_of_comparision_for_case.append(case_dict[case][0])

            mutation_matrix_dict[sv] = mutation_df

            # calculating explanation score
            ordered_col_sum = mutation_df.sum(axis=1).sort_values(ascending=False)
            ordered_row_sum = mutation_df.sum(axis=0).sort_values(ascending=False)
            explained_call_n_mosaic = ordered_col_sum[:cells_carrying_mutation_Nv].sum()
            explained_call_n_germ = ordered_row_sum[:cells_carrying_mutation_Nv].sum()
            explanation_score_mosaic = explained_call_n_mosaic / pairs_list_n
            explanation_score_germ = explained_call_n_germ / pairs_list_n
            cells_with_mutation = ",".join(list_of_cases_with_mutation)
            mutation_position = ",".join(list_of_mutation)
            comparision_for_case = ",".join(list_of_comparision_for_case)
            if cells_with_mutation == "":
                cells_with_mutation = "-"
            output_line = "\t".join(["\t".join(sv.split("\t")),
                                     str(explanation_score_mosaic),
                                     str(explanation_score_germ),
                                     str(len(list_of_cases_with_mutation)),
                                     cells_with_mutation,
                                     comparision_for_case,
                                     sv_type])
            output_file_fh.write(output_line + "\n")

            output_sv_mapping_file_fh.write(sv+"\t"+mutation_position+"\n")

        pickle.dump(mutation_matrix_dict, mutation_matrix_file_fh)
        mutation_matrix_file_fh.close()
        output_file_fh.close()
        output_sv_mapping_file_fh.close()

    def score_sv(self):
        # Extracting passed arguments
        parser = self.get_score_argument_parse_sv()
        arg = parser.parse_args(sys.argv[2:])

        # Assigning values to variable
        manifest_file = arg.manifest_file
        output_dir = arg.output_dir
        all_mutations = arg.all_mutations
        Util.ensure_dir(output_dir)

        # Extracting variant information from the manifest file.
        print("Extracting variant information")
        variant_dict, SV_mutations_dict, list_of_samples,  = self.extract_mutation_information_sv(manifest_file, output_dir,all_mutations)
        # variant_dict={mutation:[(case,control)]}
        # SV_mutations_dict={SV:{pair:[mutation]}} , mutation = chrm, start_pos, ref, alt, sv_type
        # pairs_vaf_dict={pairs:{mutation:[vaf]}}
        # list_of_samples = list of all samples in the analysis

        # Generating explanation score
        print("Generating explanation scores")
        self.explanation_score_sv(variant_dict, SV_mutations_dict, list_of_samples, output_dir)

        # Plotting
        print("Plotting")
        # plotting(output_dir,list_of_samples,af_cutoff)
        self.plot_score(output_dir)

    def plot_bar_sv(self, vaf_dict, output_dir,):
        """ plotting bar plot """
        output_file = os.path.join(output_dir, "mutation_type_count.png")
        mutation_count = {"Mosaic": 0, "Germline": 0, "Noise": 0}
        mutation_count_per_sample = {}
        variant_list = []
        for sample in vaf_dict:
            mutation_count_per_sample[sample] = {"Mosaic": 0, "Germline": 0, "Noise": 0}
            for variant_type in ["Mosaic", "Germline", "Noise"]:
                for pos, vaf, germline_score, mosaic_score in vaf_dict[sample][variant_type]:
                    if pos not in variant_list:
                        mutation_count[variant_type] += 1
                        variant_list.append(pos)
                        mutation_count_per_sample[sample][variant_type] += 1
        df = pd.DataFrame.from_dict(mutation_count, orient='index')
        title = "Number of unique mutations across all samples"
        ax = df.plot(kind='bar', legend=False, title=title, alpha=0.60, grid=False)

        for i in ax.patches:
            # get_x pulls left or right; get_height pushes up or down
            ax.text(i.get_x() + 0.15, i.get_height() / 2,
                    str(i.get_height()), fontsize=10, color='black')

        plt.ylabel("Number of mutations")
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()

        # Per sample variant counts
        df = pd.DataFrame.from_dict(mutation_count_per_sample, orient='index')
        output_file_per_sample = os.path.join(output_dir, "per_sample_mutation_count.png")
        title = "Structural variants "
        ax = df.plot(kind='bar', legend=True, grid=False)
        plt.title(title)
        plt.ylabel("Number of mutations")
        plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.5))
        plt.tight_layout()
        plt.savefig(output_file_per_sample)
        plt.close()

    def per_sample_mutation_sv(self, vaf_dict, per_sample_mutation_dir):
        """Creating per sample mutation file"""
        for sample in vaf_dict:
            file_out = os.path.join(per_sample_mutation_dir, sample+".tsv")
            file_out_fh = open(file_out,'w')
            file_out_fh.write("#SV_number\tMosaic_score\tGermline_score\tVariant_type\n")
            for variant_type in vaf_dict[sample]:
                for pos, vaf, germline_score, mosaic_score in vaf_dict[sample][variant_type]:
                    file_out_fh.write("\t".join(["\t".join(pos.split("_")), str(mosaic_score), str(germline_score), variant_type])+"\n")
            file_out_fh.close()

    def call_sv(self):
        # Extracting passed arguments
        parser = self.apply_score_argument_parse_sv()
        arg = parser.parse_args(sys.argv[2:])

        # Assigning values to variable
        get_score_dir = arg.get_score_directory
        output_dir = arg.output_dir
        mosaic_score_cutoff = float(arg.mosaic_score_cutoff)
        germline_score_cutoff = float(arg.germline_score_cutoff)

        Util.ensure_dir(output_dir)
        explanation_score_file = os.path.join(get_score_dir, "explanation_score.txt")
        vaf_dict = {}
        head = {}

        mosaic_cutoff_for_germline_mutations = float(arg.mosaic_score_cutoff_for_germline)
        germline_cutoff_for_mosaic_mutations = float(arg.germline_score_cutoff_for_mosaic)

        for i in open(explanation_score_file):
            line = i.strip().split("\t")
            if i.startswith("#"):
                for n, j in enumerate(line):
                    head[j] = n
                continue
            pos = line[head["#SV"]]
            list_of_samples = line[head["Samples_with_mutation"]].split(",")
            variant_type = ""
            germline_score = float(line[head["Germline_score"]])
            mosaic_score = float(line[head["Mosaic_score"]])

            # this block annotated the mutation as mosaic, germline or noise
            if mosaic_score > mosaic_score_cutoff:
                y_curve = y_curve = self.plot_curve([mosaic_score_cutoff], mosaic_score_cutoff, germline_score_cutoff)[0]
            else:
                y_curve = self.plot_curve([mosaic_score], mosaic_score_cutoff, germline_score_cutoff)[0]
            if germline_score > germline_score_cutoff:
                x_curve = self.plot_curve([germline_score_cutoff], germline_score_cutoff, mosaic_score_cutoff)[0]
            else:
                x_curve = self.plot_curve([germline_score], germline_score_cutoff, mosaic_score_cutoff)[0]
            if mosaic_score < x_curve and germline_score < y_curve:
                variant_type = "Noise"
            elif germline_score > y_curve and mosaic_score > x_curve and mosaic_score >= mosaic_cutoff_for_germline_mutations and germline_score > germline_cutoff_for_mosaic_mutations:
                variant_type = "Mosaic_high_freq"
            elif germline_score > y_curve and mosaic_score > x_curve and mosaic_score >= mosaic_cutoff_for_germline_mutations and germline_score <= germline_cutoff_for_mosaic_mutations:
                variant_type = "Mosaic"
            elif germline_score > y_curve and x_curve < mosaic_score < mosaic_cutoff_for_germline_mutations:
                variant_type = "Germline"
            else:
                print("variant type not detected, check for bugs in variant type definition")

            for index, sample in enumerate(list_of_samples):
                vaf = ""
                if sample in vaf_dict:
                    if variant_type in vaf_dict[sample]:
                        vaf_dict[sample][variant_type].append((pos, vaf, germline_score, mosaic_score))
                    else:
                        vaf_dict[sample][variant_type] = [(pos, vaf, germline_score, mosaic_score)]
                else:
                    vaf_dict[sample] = {}

        # structure of vaf_dict={sample:{variant_type:[(pos,vaf)]}}  ; pos = chr pos ref alt

        # per sample mutation file
        per_sample_mutation_dir = os.path.join(output_dir, "per_sample_mutation")
        Util.ensure_dir(per_sample_mutation_dir)
        self.per_sample_mutation_sv(vaf_dict,per_sample_mutation_dir)

        # plotting
        self.plot_score_annotate(explanation_score_file, output_dir, mosaic_score_cutoff, germline_score_cutoff, mosaic_cutoff_for_germline_mutations,
                                 germline_cutoff_for_mosaic_mutations)

        mutation_count_output_dir = os.path.join(output_dir, "mutation_counts")
        Util.ensure_dir(mutation_count_output_dir)
        self.plot_bar_sv(vaf_dict, mutation_count_output_dir)


if __name__ == '__main__':
    ALL2()
