using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MIC
{
    class SNP
    {
        public int id, cluster;
        public char[] genotype = new char[4500];
        public double entropy, sum, maf, x2, mi;
        /*public void init(int size){
            genotype = new char[size];
        }*/
        public SNP(int num)
        {
            Array.Resize<char>(ref genotype, num);
        }
    }

    class CENTROID
    {
        public int id, index;
        public int [] set = new int [20];
	    public double entropy, ave;
    }

    class COMBINATION
    {
        public int[] snp = new int[5];
        public double[] value = new double[4];
    }

    class BOOLEAN_EXPRESSION
    {
        public int[] snp, s_gate, i_gate, order;
        public int e_gate, TP, TN, FP, FN;
        public double mi;
        public void copy(BOOLEAN_EXPRESSION x, int num_snp)
        {
            for (int i = 0; i < num_snp; i++)
            {
                s_gate[i] = x.s_gate[i];
            }
            for (int i = 0; i < num_snp - 1; i++)
            {
                i_gate[i] = x.i_gate[i];
                order[i] = x.order[i];
            }
            e_gate = x.e_gate;
            TP = x.TP; TN = x.TN; FP = x.FP; FN = x.FN;
            mi = x.mi;
        }
        public BOOLEAN_EXPRESSION(int num_snp)
        {
            snp = new int[num_snp];
            s_gate = new int[num_snp];
            i_gate = new int[num_snp - 1];
            order = new int[num_snp - 1];
            for (int i = 0; i < num_snp - 1; i++)
            {
                s_gate[i] = 0;
                i_gate[i] = 0;
                order[i] = i;
            }
            s_gate[num_snp - 1] = 0;
        }
        public string make_string(int num_snp)
        {
            string[] boolean_expression = new string[num_snp];
            string s_temp;
            double accuracy, sensitivity, specificity, balanced_accuracy;
            for (int i = 0; i < num_snp; i++)
            {
                boolean_expression[i] = "(";
                if (0 == s_gate[i]) //and
                {
                    boolean_expression[i] += (char)('a' + i);
                    boolean_expression[i] += (char)('a' + i);
                }
                else
                {
                    boolean_expression[i] += (char)('A' + i);
                    boolean_expression[i] += (char)('a' + i);
                    boolean_expression[i] += " OR ";
                    boolean_expression[i] += (char)('a' + i);
                    boolean_expression[i] += (char)('a' + i);
                }
                boolean_expression[i] += ")";
            }
            for (int i = 0; i < num_snp - 1; i++)
            {
                if (0 == i_gate[i]) //and
                {
                    boolean_expression[order[i] + 1] = "(" + boolean_expression[i] + " AND " + boolean_expression[order[i] + 1] + ")";
                }
                else                //or
                {
                    boolean_expression[order[i] + 1] = "(" + boolean_expression[i] + " OR " + boolean_expression[order[i] + 1] + ")";
                }
            }

            if (1 == e_gate)
                boolean_expression[num_snp - 1] += "=0";
            else
                boolean_expression[num_snp - 1] += "=1";

            s_temp = "\t" + boolean_expression[num_snp - 1];
            accuracy = (double)(TP + TN) / (TP + TN + FP + FN);
            sensitivity = (double)(TP)/(TP+FN);
            specificity = (double)(TN)/(TN+FP);
            balanced_accuracy = (sensitivity+specificity)/2;
            s_temp += "\t" + accuracy;
            s_temp += "\t" + sensitivity;
            s_temp += "\t" + specificity;
            s_temp += "\t" + balanced_accuracy;
            s_temp += "\t" + TP;
            s_temp += "\t" + FP;
            s_temp += "\t" + FN;
            s_temp += "\t" + TN;

            return s_temp;
        }
        public string s_print(int num_snp)
        {
            string s_temp = "";
            /*for (int i = 0; i < num_snp; i++)
                s_temp += "\t" + s_gate[i].ToString();
            for (int i = 0; i < num_snp - 1; i++)
                s_temp += "\t" + i_gate[i].ToString() + "(" + order[i].ToString() + ")";
            if (1 == e_gate)
                s_temp += "\tnot_gate";
            s_temp += "\t" + TP.ToString() + "\t" + FN.ToString() + "\t" + FP.ToString() + "\t" + TN.ToString();*/

            s_temp = make_string(num_snp);

            return s_temp;
        }
        public int and_gate(int x, int y)
        {
            if (1 == x && 1 == y)
                return 1;
            else
                return 0;
        }
        public int or_gate(int x, int y)
        {
            if (1 == x || 1 == y)
                return 1;
            else
                return 0;
        }
        public int run(int num_snp)
        {
            for (int i = 0; i < num_snp; i++)       //each SNP
            {
                if (0 == s_gate[i])  //and
                {
                    if (0 == snp[i])
                        snp[i] = and_gate(0, 0);
                    if (1 == snp[i])
                        snp[i] = and_gate(0, 1);
                    if(2 == snp[i])
                        snp[i] = and_gate(1,1);
                }
                else        //or
                {
                    if (0 == snp[i])
                        snp[i] = or_gate(0, 0);
                    if (1 == snp[i])
                        snp[i] = or_gate(0, 1);
                    if (2 == snp[i])
                        snp[i] = or_gate(1, 1);
                }
            }
            for (int i = 0; i < num_snp - 1; i++)
            {
                if (0 == i_gate[i]) //and
                {
                    snp[order[i] + 1] = and_gate(snp[i], snp[order[i] + 1]);
                }
                else                //or
                {
                    snp[order[i] + 1] = or_gate(snp[i], snp[order[i] + 1]);
                }
            }
            return snp[num_snp-1];
        }
        public int next(int num_snp)
        {
            if (1 == next_s_gate(num_snp))
                if (1 == next_i_gate(num_snp))
                    if (1 == next_order(num_snp))
                        return 1;
            return 0;
        }
        public int next_order(int num_snp)
        {
            int index = 0;
            while (true)
            {
                if (num_snp > index + order[index]+2)
                {
                    order[index]++;
                    break;
                }
                else
                {
                    order[index] = 0;
                    index++;
                    if (num_snp-1 == index)
                        return 1;
                }
            }
            return 0;
        }
        public int next_s_gate(int num_snp)
        {
            int index = 0;
            while (true)
            {
                if (0 == s_gate[index])
                {
                    s_gate[index] = 1;
                    break;
                }
                else
                {
                    s_gate[index] = 0;
                    index++;
                    if (num_snp == index)
                        return 1;
                }
            }
            return 0;
        }
        public int next_i_gate(int num_snp)
        {
            int index = 0;
            while (true)
            {
                if (0 == i_gate[index])
                {
                    i_gate[index] = 1;
                    break;
                }
                else
                {
                    i_gate[index] = 0;
                    index++;
                    if (num_snp-1 == index)
                        return 1;
                }
            }
            return 0;
        }
    }
}
