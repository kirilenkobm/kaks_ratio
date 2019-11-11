using System;
using System.Linq;
using System.IO;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Text;

namespace kaks_ratio
{
    struct Input_data
    {
        public string first_fasta;
        public string second_fasta;
        public string first_seq_id;
        public string second_seq_id;
    }

    class Sequences
    {
        private readonly char[] letters = { 'A', 'T', 'G', 'C' };
        public readonly string[] stops = { "TAA", "TAG", "TGA" };

        public string Extract_sequence
            (string fasta_file, string seq_id)
        {
            string file_content = System.IO.File.ReadAllText(fasta_file);
            string[] splitted = file_content.Split(">");
            string seq = "";

            foreach (var elem in splitted)
            {
                string[] elem_split = elem.Split("\n");
                string this_seq_id = elem_split[0];
                if (!this_seq_id.Equals(seq_id)) { continue; }
                // this is equal
                string[] seq_elems = elem_split.Skip(1).ToArray();
                seq = string.Join("", seq_elems);
                return seq;
            }
            // if we are here -> the sequence we need was not found
            Console.WriteLine($"Error! Sequence {seq_id} not found in {fasta_file}");
            Console.WriteLine("Abort.");
            Environment.Exit(1);
            return seq;
        }

        public IEnumerable<string> Parts(string str, int chunkSize)
        {
            return Enumerable.Range(0, str.Length / chunkSize)
                .Select(i => str.Substring(i * chunkSize, chunkSize));
        }

        public string get_AA(string codon)
        {
            switch (codon)
            {
                case "ATG":
                    return "M";
                case "ATT":
                case "ATC":
                case "ATA":
                    return "I";
                case "GTT":
                case "GTC":
                case "GTA":
                case "GTG":
                    return "V";
                case "CTT":
                case "CTC":
                case "CTA":
                case "CTG":
                case "TTA":
                case "TTG":
                    return "L";
                case "TTT":
                case "TTC":
                    return "F";
                case "TCT":
                case "TCC":
                case "TCA":
                case "TCG":
                    return "S";
                case "CCT":
                case "CCC":
                case "CCA":
                case "CCG":
                    return "P";
                case "ACT":
                case "ACC":
                case "ACA":
                case "ACG":
                    return "T";
                case "GCT":
                case "GCC":
                case "GCA":
                case "GCG":
                    return "A";
                case "TAT":
                case "TAC":
                    return "Y";
                case "TGT":
                case "TGC":
                    return "C";
                case "CAU":
                case "CAC":
                    return "H";
                case "CAA":
                case "CAG":
                    return "Q";
                case "CGT":
                case "CGC":
                case "CGA":
                case "CGG":
                case "AGA":
                case "AGG":
                    return "R";
                case "AGT":
                case "AGC":
                    return "S";
                case "AAT":
                case "AAC":
                    return "N";
                case "AAA":
                case "AAG":
                    return "K";
                case "GAT":
                case "GAC":
                    return "D";
                case "GAA":
                case "GAG":
                    return "E";
                case "GGT":
                case "GGC":
                case "GGA":
                case "GGG":
                    return "G";
                case "TAA":
                case "TAG":
                case "TGA":
                    return "*";
                case "TGG":
                    return "W";
                default:
                    return "X";
            }
        }

        private Tuple<char, char, char> codon_to_tuple(string codon)
        {
            char num_1 = codon[0];
            char num_2 = codon[1];
            char num_3 = codon[2];
            return Tuple.Create(num_1, num_2, num_3);
        }


        private string tuple_to_codon(Tuple<char, char, char> codon_tup)
        {
            string ans = "";
            ans.Append(codon_tup.Item1);
            ans.Append(codon_tup.Item2);
            ans.Append(codon_tup.Item3);
            return ans;
        }


        public Tuple<UInt32, UInt32> codon_sites_count(string codon)
        {
            if (codon.Contains("N"))
            {
                // we cannot count anything
                return Tuple.Create((UInt32)0, (UInt32)0);
            }
            // no N, we can count variants
            UInt32 syn_sites = 0;
            UInt32 nsyn_sites = 0;
            string codon_AA = get_AA(codon);

            for (int i = 0; i < 3; ++i)
            {
                char c = codon[i];
                char[] not_c = letters.Where(w => w != c).ToArray();
                foreach (var nc in not_c)
                {
                    StringBuilder var_codon = new StringBuilder(codon);
                    var_codon[i] = nc;
                    string var_AA = get_AA(var_codon.ToString());
                    if (codon_AA.Equals(var_AA))
                    {
                        ++syn_sites;
                    } else
                    {
                        ++nsyn_sites;
                    }
                }
            }
    
            return Tuple.Create(syn_sites, nsyn_sites);
        }

        public UInt32 get_codon_dist(string codon_1, string codon_2)
        {
            UInt32 dist = 0;
            for (int i = 0; i < 3; ++i)
            {
                if (codon_1[i] != codon_2[i])
                {
                    ++dist;
                }
            }
            return dist;
        }
    }


    class kaks_ratio
    {

        static void Main(string[] args)
        {
            if (args.Length < 4
                || args.Contains("--help")
                || args.Contains("-h"))
            {
                Console.Write($"Usage: {AppDomain.CurrentDomain.FriendlyName}");
                Console.Write(" [first fasta] [first sequence name]");
                Console.WriteLine(" [second fasta] [second sequence name]");
                Console.WriteLine("If second fasta is the same with the first, just use - ");
                Console.WriteLine("The same rule works for the second sequence id.");
                Environment.Exit(0);
            }

            // write input data to struct
            Input_data input_data = new Input_data();
            input_data.first_fasta = args[0];
            input_data.first_seq_id = args[1];

            // second elems are bit more tricky
            if (args[2].Equals("-"))
            {
                input_data.second_fasta = args[0];
            } else
            {
                input_data.second_fasta = args[2];
            }
            if (args[3].Equals("-"))
            {
                input_data.second_seq_id = args[1];
            } else
            {
                input_data.second_seq_id = args[3];
            }

            // read first and second sequences
            Sequences seq_tools = new Sequences();
            string first_seq = seq_tools.Extract_sequence(input_data.first_fasta,
                                                          input_data.first_seq_id);
            string second_seq = seq_tools.Extract_sequence(input_data.second_fasta,
                                                           input_data.second_seq_id);
            if (first_seq.Length != second_seq.Length)
            {
                Console.WriteLine("Error! Sequences must have the same length!");
                Console.WriteLine("Otherwise the computation makes no sense.");
                Console.WriteLine("Abort.");
                Environment.Exit(1);
            } else if (first_seq.Length % 3 != 0)
            {
                Console.WriteLine("Error! Sequences length mod 3 must be 0!");
                Console.WriteLine("A codon alignment required for the computations.");
                Console.WriteLine("Abort.");
                Environment.Exit(1);
            }

            string[] first_codons = seq_tools.Parts(first_seq, 3).ToArray();
            string[] second_codons = seq_tools.Parts(second_seq, 3).ToArray();

            UInt32 diff_p_counter = 1;
            UInt32 seq_changes = 0;
            UInt32 codon_changes = 0;
            UInt32 total_syn_codons = 0;
            UInt32 total_nsyn_codons = 0;

            UInt16 syn_muts = 0;
            UInt16 non_syn_muts = 0;

            for (int i = 0; i < first_codons.Length; ++i)
            {
                string codon_1 = first_codons[i];
                string codon_2 = second_codons[i];

                // if so -> cannot count anything
                if (codon_1.Contains("N") || codon_2.Contains("N")) { continue; }
                // if stops -> also no way to compute anything
                if (seq_tools.stops.Contains(codon_1)
                    || seq_tools.stops.Contains(codon_2)) { continue; }

                Tuple<UInt32, UInt32> codon_1_sites = seq_tools.codon_sites_count(codon_1);
                total_syn_codons += codon_1_sites.Item1;
                total_nsyn_codons += codon_1_sites.Item2;

                if (codon_1.Equals(codon_2)) {continue;}
                // if we're here -> codons differ

                string AA_1 = seq_tools.get_AA(codon_1);
                string AA_2 = seq_tools.get_AA(codon_2);
                UInt32 codon_diff = seq_tools.get_codon_dist(codon_1, codon_2);
                seq_changes += codon_diff;
                ++codon_changes;
                string change_descr;

                if (AA_1.Equals(AA_2))
                {
                    ++syn_muts;
                    change_descr = "SYN";
                } else
                {
                    ++non_syn_muts;
                    change_descr = "NON_SYN";
                }

                Console.Write($"{diff_p_counter} | CODON NUM {i + 1} | {codon_1} -> {codon_2} ");
                Console.Write($"{AA_1} -> {AA_2} {change_descr}\n");
                ++diff_p_counter;
            }

            // compute omega
            double omega_base;
            double non_syn_to_syn;
            double omega;

            if (total_syn_codons > 0)
            {
                // it syn codons is 0 -> zerodivision error shall be
                omega_base = (double)total_nsyn_codons / (double)total_syn_codons;
            } else { omega_base = 99999.9; }

            if (syn_muts > 0) { non_syn_to_syn = (double)non_syn_muts / (double)syn_muts; }
            else { non_syn_to_syn = 99999.9; }

            if (syn_muts > 0)
            {
                omega = non_syn_to_syn / omega_base;
            }else { omega = 99999.9; }

            // write the report
            Console.WriteLine($"Overall {codon_changes} different codons and {seq_changes} changes.");
            Console.WriteLine($"{syn_muts} synonymous and {non_syn_muts} non-synonymous changes");
            Console.WriteLine($"There are {total_syn_codons} synonymous sites and {total_nsyn_codons} non-synonymous.");
            Console.WriteLine($"Syn/non_syn sites ratio is {omega_base}");
            Console.WriteLine($"Omega = {omega}");
        }
    }
}
