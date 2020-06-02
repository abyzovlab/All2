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


class ALL2():
    def __init__(self):
        parser = argparse.ArgumentParser(description='All to all comparision',
                                         usage=''' python ALL2.py <command> [<args>]
        Three commands to use:
                score --> Generates mosaic and germline scores.
                call --> Based on score cut-off generates sample level files/plots for mosaic,
                        germline mutations and plots variant allele frequency, mutation spectrum
                matrix --> Plot the mutation matrix
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
        parser.add_argument("-m", "--manifest_file", help="Path to manifest file", required=True,
                            type=Util.FileValidator)
        parser.add_argument("-o", "--output_dir", help="Path to directory where results will be written",
                            required=True)
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
        parser.add_argument("-ms", "--mosaic_score_cutoff", help="Mosaic score cut-off (default=1)", default="1.0")
        parser.add_argument("-gs", "--germline_score_cutoff", help="Germline score cut-off (default=1)", default="1.0")
        return parser

    def matrix(self):
        parser = self.mutation_matrix_plot_argument_parse()
        arg = parser.parse_args(sys.argv[2:])
        import seaborn as sns
        # Assigning values to variable
        ALL2_output = arg.get_score_directory
        output_dir = arg.output_dir
        mutation_list = arg.mutation
        explaination_file = os.path.join(ALL2_output, "explaination_score.txt")
        Util.ensure_dir(output_dir)

        explaination_dict = {}
        head = {}
        for i in open(explaination_file):
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
            mutation = "_".join([chrm, pos, ref, alt])
            mutation_related_info = {"mosaic_score": mosaic_score, "germline_score": germline_score,
                                     "sample": samples, "vaf_samples": vaf_samples}
            explaination_dict[mutation] = mutation_related_info
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

            list_of_samples = explaination_dict[mutation]["sample"]
            list_of_samples_vaf = explaination_dict[mutation]["vaf_samples"]
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
            ax1.set_title(mutation)
            ax2.set_title("VAF", fontsize=10)
            # mosaic_score = explaination_dict[mutation]["mosaic_score"]
            # germline_score = explaination_dict[mutation]["germline_score"]
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, mutation + ".png"))
            plt.close()

    def extract_mutation_information(self, manifest_file, output_dir):
        variant_dict = {}
        pairs_vaf_dict = {}
        head = {}
        for i in open(manifest_file):
            line = i.strip().split("\t")
            if i.startswith("#"):
                for n, j in enumerate(line):
                    head[j.replace("#", "")] = n
                continue
            case = line[head["Case"]]
            control = line[head["Control"]]
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
                variant = variant.decode("utf-8")
                line = variant.strip().split("\t")
                if variant.startswith("#"):
                    if variant.startswith("#CHROM"):
                        for n, j in enumerate(line):
                            variant_head[j] = n
                        if line[variant_head["FORMAT"] + 1].upper() == "TUMOR":
                            case = line[variant_head["FORMAT"] + 1]
                        elif line[variant_head["FORMAT"] + 2].upper() == "TUMOR":
                            case = line[variant_head["FORMAT"] + 2]
                    continue
                chrm = line[variant_head["#CHROM"]]
                pos = line[variant_head["POS"]]
                ref = line[variant_head["REF"]]
                all_alt = line[variant_head["ALT"]]
                # Getting AD and DP field for case
                case_format = line[variant_head["FORMAT"]].split(":")
                try:
                    case_genotype = line[variant_head[case]].split(":")
                except KeyError:
                    print("Please make sure the name of the case and control in the manifest file match the"
                          " case and control specified in the vcf")
                    exit()
                try:
                    ad_index = case_format.index("AD")
                except ValueError:
                    case_genotype_ad = "Absent"
                try:
                    depth_index = case_format.index("DP")
                except ValueError:
                    case_genotype_depth = "Absent"

                if depth_index is not "Absent":
                    case_genotype_depth = case_genotype[depth_index]
                if ad_index is not "Absent":
                    case_genotype_ad = case_genotype[ad_index]
                if case_genotype_depth == "0":
                    continue
                # making sure to take care of multiallelic locations
                for n, alt in enumerate(all_alt.split(",")):
                    mutation = "_".join([chrm, pos, ref, alt])
                    # storing variants
                    if mutation in variant_dict:
                        variant_dict[mutation].append(pair)
                    else:
                        variant_dict[mutation] = [pair]
                    # storing VAFs
                    if case_genotype_ad is not "Absent":
                        if case_genotype_depth is not "Absent":
                            alt_supporting_read = case_genotype_ad.split(",")[n + 1]
                            vaf = str(float(float(alt_supporting_read) / float(case_genotype_depth)))
                        else:
                            vaf = "Absent"
                    else:
                        vaf = "Absent"
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
        return variant_dict, pairs_vaf_dict, list_of_samples

    def explaination_score(self, variant_dict, pairs_vaf_dict, list_of_samples, output_dir):

        output_file = os.path.join(output_dir, "explaination_score.txt")
        output_file_fh = open(output_file, 'w')
        output_file_fh.write("#Chrom\tPos\tRef\tAlt\tMosaic_score\tGermline_score\tNumber_of_samples_with_mutation"
                             "\tSamples_with_mutation\tVAF_of_samples_with_mutation\tNumber_of_comparision_per_sample"
                             "\n")

        # number_of_cells_N is the total number of cells in the experiment
        total_number_of_cells_N = len(list_of_samples)

        # Preparing to store the data matrix for each mutation
        mutation_matrix_dict = {}
        mutation_matrix_file = os.path.join(output_dir, "mutation_matrix.pkl")
        mutation_matrix_file_fh = open(mutation_matrix_file, 'wb')

        for mutation in variant_dict:
            # pairs_list is a list ofpairs the mutation was called in
            pairs_list = variant_dict[mutation]
            # pairs_list_n is the number of pairs the mutation was called in
            pairs_list_n = len(pairs_list)
            # cell_fraction_f is the fraction of cells carrying the mutatin
            cell_fraction_f = 1 / 2 - sqrt(1 / 4 - pairs_list_n / total_number_of_cells_N ** 2)
            # is the number of cells carrying the mutation
            cells_carrying_mutation_Nv = round(cell_fraction_f * total_number_of_cells_N)
            # Creating an 'zero' data frame/matrix and updating mutation specific dataframe/matrix
            mutation_df = pd.DataFrame(np.zeros((total_number_of_cells_N, total_number_of_cells_N)),
                                       index=list_of_samples, columns=list_of_samples)
            list_of_cases_with_mutation = []
            list_of_vaf_cases_with_mutation = []
            list_of_comparision_for_case = []
            case_dict = {}
            for case, control in pairs_list:
                vaf = str(pairs_vaf_dict[(case, control)][mutation][0])
                if case in case_dict:
                    case_dict[case][0] = vaf
                    case_dict[case][1] = str(int(case_dict[case][1]) + 1)
                else:
                    case_dict[case] = [vaf, "1"]

                mutation_df.loc[case, control] = 1

            for case in case_dict:
                list_of_cases_with_mutation.append(case)
                list_of_vaf_cases_with_mutation.append(case_dict[case][0])
                list_of_comparision_for_case.append(case_dict[case][1])

            mutation_matrix_dict[mutation] = mutation_df

            # calculating explaination score
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
            output_line = "\t".join(["\t".join(mutation.split("_")),
                                     str(explanation_score_mosaic),
                                     str(explanation_score_germ),
                                     str(len(list_of_cases_with_mutation)),
                                     cells_with_mutation,
                                     vaf_of_samples_with_mutation,
                                     comparision_for_case])
            output_file_fh.write(output_line + "\n")

        pickle.dump(mutation_matrix_dict, mutation_matrix_file_fh)
        mutation_matrix_file_fh.close()
        output_file_fh.close()

    def plot_score(self, output_dir):
        explaination_score = os.path.join(output_dir, "explaination_score.txt")

        # Plotting germline versus mosaic score scatter plot
        df_es = pd.read_table(explaination_score)
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
        df_es.plot.scatter("Mosaic_score", "Germline_score", c=size, s=[float(s / 10.0) for s in size], cmap='tab10')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "Explaination_score_scatter.png"))
        plt.close()

    def score(self):
        # Extracting passed arguments
        parser = self.get_score_argument_parse()
        arg = parser.parse_args(sys.argv[2:])

        # Assigning values to variable
        manifest_file = arg.manifest_file
        output_dir = arg.output_dir
        Util.ensure_dir(output_dir)

        # Extracting variant information from the manifest file.
        print("Extracting variant information")
        variant_dict, pairs_vaf_dict, list_of_samples = self.extract_mutation_information(manifest_file, output_dir)
        # variant_dict={mutation:[(case,control)]}
        # pairs_vaf_dict={pairs:{mutation:[vaf]}}
        # list_of_samples = list of all samples in the analysis

        # Generating explaination score
        print("Generating explaination scores")
        self.explaination_score(variant_dict, pairs_vaf_dict, list_of_samples, output_dir)

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
                for pos, vaf in vaf_dict[sample][variant_type]:
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
                for pos, vaf in vaf_dict[sample][variant_type]:
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
                                      alpha=0.5, legend=None,
                                      weights=np.ones_like(df_snp[df_snp.columns[0]]) * 1. / len(df_snp),
                                      title=title)
                    df_indel.plot(kind='hist', bins=33, range=(0, 1),
                                  alpha=0.5, legend=None,
                                  weights=np.ones_like(df_indel[df_indel.columns[0]]) * 1. / len(df_indel),
                                  ax=ax1)
                    plt.legend(["SNP", "INDEL"])
                elif df_indel.empty:
                    df_snp.plot(kind='hist', bins=33, range=(0, 1),
                                alpha=0.5, legend=None,
                                weights=np.ones_like(df_snp[df_snp.columns[0]]) * 1. / len(df_snp),
                                title=title)
                    plt.legend(["SNP"])
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
        for pos, vaf in vaf_dict[sample][variant_type]:
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
        explaination_score_file = os.path.join(get_score_dir, "explaination_score.txt")
        vaf_dict = {}
        head = {}
        for i in open(explaination_score_file):
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
            if mosaic_score >= mosaic_score_cutoff:
                variant_type = "Mosaic"
            elif germline_score >= germline_score_cutoff:
                variant_type = "Germline"
            else:
                variant_type = "Noise"
            for index, sample in enumerate(list_of_samples):
                vaf = float(list_of_vafs[index])
                if sample in vaf_dict:
                    if variant_type in vaf_dict[sample]:
                        vaf_dict[sample][variant_type].append((pos, vaf))
                    else:
                        vaf_dict[sample][variant_type] = [(pos, vaf)]
                else:
                    vaf_dict[sample] = {}

        # structure of var_dict={sample:{variant_type:[(pos,vaf)]}}

        # plotting
        mutation_count_output_dir = os.path.join(output_dir, "mutation_counts")
        Util.ensure_dir(mutation_count_output_dir)
        self.plot_bar(vaf_dict, af_cutoff, mutation_count_output_dir)

        vaf_output_dir = os.path.join(output_dir, "vaf_plots")
        Util.ensure_dir(vaf_output_dir)
        self.plot_vaf(vaf_dict, af_cutoff, vaf_output_dir)

        mutation_spectrum_output_dir = os.path.join(output_dir, "mutation_spectrum")
        Util.ensure_dir(mutation_spectrum_output_dir)
        self.plot_mutation_spectrum(vaf_dict, af_cutoff, mutation_spectrum_output_dir, reference)


if __name__ == '__main__':
    ALL2()
