using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Collections;
using System.Threading;
//using Meta.Numerics;
using MathNet.Numerics.Distributions;
using Microsoft.Office.Interop.Excel;

namespace MIC
{
    public partial class Form1 : Form
    {
        const double TO_ZERO = 0.000000001;
        const double LN2 = 0.693147181;
        int num_case;
        int num_control;
        int num_sample;
        SNP[] snp = new SNP[1024]; 
        int num_SNP = 0, set_size = 10;
        double h_type;
        CENTROID [] centroid = new CENTROID[2];

        public Form1()
        {
            InitializeComponent();
        }

        //case file selection
        private void button1_Click(object sender, EventArgs e)
        {
            openFileDialog1.ShowDialog();
        }   

        private void button2_Click(object sender, EventArgs e)
        {
            openFileDialog2.ShowDialog();
        }

        private void openFileDialog1_FileOk(object sender, CancelEventArgs e)
        {
            textBox1.Text = openFileDialog1.FileName;
        }

        private void openFileDialog2_FileOk(object sender, CancelEventArgs e)
        {
            textBox2.Text = openFileDialog2.FileName;
        }

        private void loadData()
        {
            CheckForIllegalCrossThreadCalls = false;
            label9.Text = "Data loading";
            label9.Refresh();
            try
            {
                FileStream read_file1 = new FileStream(textBox1.Text, FileMode.Open);
                StreamReader inStream1 = new StreamReader(read_file1, System.Text.Encoding.Default);
                FileStream read_file2 = new FileStream(textBox2.Text, FileMode.Open);
                StreamReader inStream2 = new StreamReader(read_file2, System.Text.Encoding.Default);
                string temp = "Init";
                int i = 0;
                while(true)
                //for (int i = 0; i < num_SNP_orin; i++)
                {
                    temp = inStream1.ReadLine();
                    if (null == temp)
                        break;
                    num_case = temp.Length;     //added
                    char[] data_case = new char[num_case];
                    //int[] id_case = new int[num_SNP_orin];
                    //string[] temp2 = temp.Split(new char[] { ' ', '\t', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    //id_case[i] = i + 1; // Convert.ToInt32(temp2[0]);
                    //temp = temp2[1];
                    
                    for (int j = 0; j < num_case; j++)
                        data_case[j] = (char)temp[j];
                    temp = inStream2.ReadLine();
                    num_control = temp.Length;
                    char[] data_control = new char[num_control];
                    //int[] id_control = new int[num_SNP_orin];

                    num_sample = num_case + num_control;

                    h_type = -(double)num_case / num_sample * Math.Log((double)num_case / num_sample) / Math.Log(2.0)
                        - (double)num_control / num_sample * Math.Log((double)num_control / num_sample) / Math.Log(2.0);

                    //string[] temp3 = temp.Split(new char[] { ' ', '\t', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    //id_control[i] = i + 1;  // Convert.ToInt32(temp2[0]);
                    //temp = temp3[1];
                                                           
                    for (int j = 0; j < num_control; j++)
                        data_control[j] = (char)temp[j];
                    
                    //cal_entropy
                    int cnt = 0, cnt_case = 0, cnt_control = 0;
                    int[] table = { 0, 0, 0 };
                    double[] p_index = { 0, 0, 0 };
                    double[,] p_index1 = { { 0, 0, 0 }, { 0, 0, 0 } };
                    double maf, missing;
                    for (int j = 0; j < num_case; j++)
                    {
                        if ('0' <= data_case[j] && data_case[j] <= '2')
                        {
                            table[data_case[j] - '0']++;
                            p_index1[0, data_case[j] - '0']++;
                            cnt++;
                            cnt_case++;
                        }
                    }
                    for (int j = 0; j < num_control; j++)
                    {
                        if ('0' <= data_control[j] && data_control[j] <= '2')
                        {
                            table[data_control[j] - '0']++;
                            p_index1[1, data_control[j] - '0']++;
                            cnt++;
                            cnt_control++;
                        }
                    }
                    for (int j = 0; j < 3; j++)
                    {
                        if (table[j] > 0)
                            p_index[j] = (double)table[j] / cnt;
                        else
                            p_index[j] = TO_ZERO;

                        if (p_index1[0, j] > 0)
                            p_index1[0, j] /= cnt;
                        else
                            p_index1[0, j] = TO_ZERO;

                        if (p_index1[1, j] > 0)
                            p_index1[1, j] /= cnt;
                        else
                            p_index1[1, j] = TO_ZERO;
                    }
                    if (p_index[0] > p_index[2])
                        maf = 2 * p_index[2] + p_index[1];
                    else
                        maf = 2 * p_index[0] + p_index[1];
                    maf /= 2;
                    missing = (double)(num_sample - cnt) / num_sample;
                    if (maf > 0.05 && missing < 0.05)
                    {
                        snp[num_SNP] = new SNP(num_case+num_control);
                        snp[num_SNP].entropy = 0;
                        snp[num_SNP].mi = 0;
                        snp[num_SNP].id = i+1;
                        
                        for (int j = 0; j < num_case; j++)
                            snp[num_SNP].genotype[j] = data_case[j];
                        for (int j = 0; j < num_control; j++)
                            snp[num_SNP].genotype[num_case + j] = data_control[j];
                        for (int j = 0; j < 3; j++)
                        {
                            snp[num_SNP].entropy += (-(p_index[j]) * (Math.Log(p_index[j]) / Math.Log(2.0)));
                            snp[num_SNP].mi += (p_index1[0, j] * (Math.Log(p_index1[0, j]) / Math.Log(2.0)));
                            snp[num_SNP].mi += (p_index1[1, j] * (Math.Log(p_index1[1, j]) / Math.Log(2.0)));
                        }
                        snp[num_SNP].mi += (snp[num_SNP].entropy + h_type);
                        snp[num_SNP].maf = maf;
                        num_SNP++;
                        if (0 == num_SNP % 1000)
                        {
                            progressBar1.Value = (num_SNP / 1000)%101;
                            label9.Text = num_SNP.ToString() + " SNPs";
                            label9.Refresh();
                        }
                        if (0 == num_SNP % 1024)
                            Array.Resize<SNP>(ref snp, num_SNP + 1024);
                    }
                    i++;
                }

                textBox5.Text = num_case.ToString();
                textBox5.Refresh();
                textBox7.Text = num_control.ToString();
                textBox7.Refresh();
                textBox4.Text = i.ToString();
                textBox4.Refresh();

                Array.Resize<SNP>(ref snp, num_SNP);
                inStream1.Close();
                read_file1.Close();
                inStream2.Close();
                read_file2.Close();
            }
            catch (IOException e_IO)
            {
                label9.Text = "Loading error" + e_IO.Message;
            }
            label9.Text = num_SNP + " SNPs are remianed";
            label9.Refresh();
        }             
        private static int snpCompare(SNP x, SNP y)
        {
            if (x.mi < y.mi)
                return 1;
            else
                return -1;
        }
        private static int combiCompare_0(COMBINATION x, COMBINATION y)
        {
            if (x.value[0] < y.value[0])
                return 1;
            else
                return -1;
        }
        private static int combiCompare_1(COMBINATION x, COMBINATION y)
        {
            if (x.value[1] < y.value[1])
                return 1;
            else
                return -1;
        }
        private static int combiCompare_2(COMBINATION x, COMBINATION y)
        {
            if (x.value[2] < y.value[2])
                return 1;
            else
                return -1;
        }
        private static int combiCompare_3(COMBINATION x, COMBINATION y)
        {
            if (x.value[3] < y.value[3])
                return 1;
            else
                return -1;
        }

        private double distance_2(int a, int b, int opt)
        {
            double h_xy = 0, h_xyt=0, h_xt=0, h_yt=0;
	        double[,] table = new double [3,3];
            double[,] t_case = new double [3,3];
            double[,] t_control = new double [3,3];
            double [] px = new double[3], py = new double[3];
            double[] px_case = new double[3], px_control = new double[3], py_case = new double[3], py_control = new double[3];
	        double P_AB, P_A=0, P_B=0, P_a=0, P_b=0;
	        double t_ratio;
	        int cnt = 0;
	        double x2;
	        for(int i=0; i<3;i++){
		        for(int j=0;j<3;j++){
			        table[i,j] = 0;
			        t_case[i,j] = 0;
			        t_control[i,j] = 0;
		        }
		        px[i]=0;
		        py[i]=0;
                px_case[i] = 0;
                px_control[i] = 0;
                py_case[i] = 0;
                py_control[i] = 0;
	        }	        

	        for(int i=0; i<num_sample; i++)
	        {
		        if(snp[a].genotype[i]-'0' >= 0 && snp[b].genotype[i] -'0' >=0)
		        {
			        table[snp[a].genotype[i]-'0', snp[b].genotype[i]-'0']++;
			        
			        cnt++;
			        if(i<num_case)
			            t_case[snp[a].genotype[i]-'0', snp[b].genotype[i]-'0']++;
				    else
			            t_control[snp[a].genotype[i]-'0', snp[b].genotype[i]-'0']++;
		        }
	        }

	        t_ratio = (double)num_case/num_control;
	        for(int i=0; i<3; i++)
		        for(int j=0;j<3; j++)
		        {
			        table[i,j] /= cnt;
			        if(table[i,j] > TO_ZERO)
				        h_xy += -table[i,j] * Math.Log(table[i,j])/Math.Log(2.0);

			        t_case[i,j] /= cnt;
			        t_control[i,j] /= cnt;
			        if(t_case[i,j] > TO_ZERO)
				        h_xyt += -t_case[i,j] * Math.Log(t_case[i,j])/Math.Log(2.0);
			        if(t_control[i,j] > TO_ZERO)
				        h_xyt += -t_control[i,j] * Math.Log(t_control[i,j])/Math.Log(2.0);						
		        }

	        P_AB = table[0,0] + (table[0,1] + table[1,0])/2 + table[1,1]/4;
	
	        x2 = 0;
	        for(int i=0;i<3; i++){
		        P_A += table[0,i] + table[1,i]/2;
		        P_a += table[2,i] + table[1,i]/2;
		        P_B += table[i,0] + table[i,1]/2;
		        P_b += table[i,2] + table[i,1]/2;
                for(int j=0;j<3;j++)
		        {
                    px_case[i] += t_case[i, j];
                    px_control[i] += t_control[i, j];
                    py_case[i] += t_case[j, i];
                    py_control[i] += t_control[j, i];

			        px[i] += t_case[i,j] + t_control[i,j];
			        py[j] += t_case[i,j] + t_control[i,j];			       
		        }
                
		        if(px_case[i] > TO_ZERO)
			        h_xt += -px_case[i]*Math.Log(px_case[i])/Math.Log(2.0);
		        if(px_control[i] > TO_ZERO)
			        h_xt += -px_control[i]*Math.Log(px_control[i])/Math.Log(2.0);
		        if(py_case[i] > TO_ZERO)
			        h_yt += -py_case[i]*Math.Log(py_case[i])/Math.Log(2.0);
		        if(py_control[i] > TO_ZERO)
			        h_yt += -py_control[i]*Math.Log(py_control[i])/Math.Log(2.0);
	        }
	        for(int i=0;i<3;i++)
                for (int j = 0; j < 3; j++)
                {
                    double e_case, e_control; //for x2
                    e_case = (px[i] * py[j] * num_case + TO_ZERO);
                    e_control = (px[i] * py[j] * num_control + TO_ZERO);
                    t_case[i, j] *= num_sample;
                    t_control[i, j] *= num_sample;
                    x2 += (t_case[i, j] - e_case) * (t_case[i, j] - e_case) / (e_case + TO_ZERO);
                    x2 += (t_control[i, j] - e_control) * (t_control[i, j] - e_control) / (e_control + TO_ZERO);
                }
	
	        if (1 == opt) //X2
		        return x2;	//임상섭
            else
                return (-h_xyt+h_type+h_xy);	//임상섭
        }
        private double distance_3(int a, int b, int c, int opt)
        {
            double h_xy = 0, h_xyt = 0;
            double[,,] table = new double[3, 3, 3];
            double[,,] t_case = new double[3, 3, 3];
            double[,,] t_control = new double[3, 3, 3];
            double[] pa_case = new double[3], pa_control = new double[3], pb_case = new double[3], pb_control = new double[3], pc_case = new double[3], pc_control = new double[3];
           
            double t_ratio;
            int cnt = 0, cnt_case=0, cnt_control = 0;
            double x2;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        table[i, j, k] = 0;
                        t_case[i, j, k] = 0;
                        t_control[i, j, k] = 0;
                    }
                }
                pa_case[i] = 0;
                pa_control[i] = 0;
                pb_case[i] = 0;
                pb_control[i] = 0;
                pc_case[i] = 0;
                pc_control[i] = 0;
            }

            for (int i = 0; i < num_sample; i++)
            {
                if (snp[a].genotype[i] - '0' >= 0 && snp[b].genotype[i] - '0' >= 0 && snp[c].genotype[i] - '0' >= 0)
                {
                    table[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0']++;
                    cnt++;
                    if (i < num_case)
                    {
                        t_case[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0']++;
                        cnt_case++;
                    }
                    else
                    {
                        t_control[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0']++;
                        cnt_control++;
                    }
                }
            }

            t_ratio = (double)num_case / num_control;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    for(int k=0; k<3; k++)
                    {
                        table[i, j, k] /= cnt;
                        if (table[i, j, k] > TO_ZERO)
                            h_xy += -table[i, j, k] * Math.Log(table[i, j, k]) / Math.Log(2.0);

                        t_case[i, j, k] /= cnt;
                        t_control[i, j, k] /= cnt;
                        if (t_case[i, j, k] > TO_ZERO)
                            h_xyt += -t_case[i, j, k] * Math.Log(t_case[i, j, k]) / Math.Log(2.0);
                        if (t_control[i, j, k] > TO_ZERO)
                            h_xyt += -t_control[i, j, k] * Math.Log(t_control[i, j, k]) / Math.Log(2.0);
                    }
                }

            x2 = 0;
	        for(int i=0;i<3; i++){
		        for(int j=0;j<3;j++)
			        for(int k=0;k<3;k++)
			        {				
				        double e_case, e_control; //for x2
				        e_case=(t_case[i,j,k] + t_control[i,j,k] + TO_ZERO)*cnt_case;
				        e_control=(t_case[i,j,k] + t_control[i,j,k] + TO_ZERO)*cnt_control;
				        x2 += (t_case[i,j,k] * cnt - e_case)*(t_case[i,j,k] * cnt - e_case)/(e_case+TO_ZERO);
				        x2 += (t_control[i,j,k] * cnt-e_control) * (t_control[i,j,k] *cnt-e_control)/(e_control+TO_ZERO);
			        }
	        }

            if (1 == opt) //X2
                return x2;	//임상섭
            else
                return (-h_xyt + h_type + h_xy);	//임상섭
        }
        private double distance_4(int a, int b, int c, int d, int opt)
        {
            double h_xy = 0, h_xyt = 0;
            double[, ,,] table = new double[3, 3, 3,3];
            double[, ,,] t_case = new double[3, 3, 3,3];
            double[, ,,] t_control = new double[3, 3, 3,3];
            double[] pa_case = new double[3], pa_control = new double[3], 
                pb_case = new double[3], pb_control = new double[3],
                pc_case = new double[3], pc_control = new double[3],
                pd_case = new double[3], pd_control = new double[3];

            double t_ratio;
            int cnt = 0, cnt_case = 0, cnt_control = 0;
            double x2;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            table[i, j, k, l] = 0;
                            t_case[i, j, k, l] = 0;
                            t_control[i, j, k, l] = 0;
                        }
                    
                pa_case[i] = 0;
                pa_control[i] = 0;
                pb_case[i] = 0;
                pb_control[i] = 0;
                pc_case[i] = 0;
                pc_control[i] = 0;
                pd_case[i] = 0;
                pd_control[i] = 0;
            }

            for (int i = 0; i < num_sample; i++)
            {
                if (snp[a].genotype[i] - '0' >= 0 && snp[b].genotype[i] - '0' >= 0 && snp[c].genotype[i] - '0' >= 0 && snp[d].genotype[i] - '0' >= 0)
                {
                    table[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0', snp[d].genotype[i] - '0']++;
                    cnt++;
                    if (i < num_case)
                    {
                        t_case[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0', snp[d].genotype[i] - '0']++;
                        cnt_case++;
                    }
                    else
                    {
                        t_control[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0', snp[d].genotype[i] - '0']++;
                        cnt_control++;
                    }
                }
            }

            t_ratio = (double)num_case / num_control;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            table[i, j, k,l] /= cnt;
                            if (table[i, j, k, l] > TO_ZERO)
                                h_xy += -table[i, j, k, l] * Math.Log(table[i, j, k, l]) / Math.Log(2.0);

                            t_case[i, j, k, l] /= cnt;
                            t_control[i, j, k, l] /= cnt;
                            if (t_case[i, j, k, l] > TO_ZERO)
                                h_xyt += -t_case[i, j, k, l] * Math.Log(t_case[i, j, k, l]) / Math.Log(2.0);
                            if (t_control[i, j, k, l] > TO_ZERO)
                                h_xyt += -t_control[i, j, k, l] * Math.Log(t_control[i, j, k, l]) / Math.Log(2.0);
                        }                        
            

            x2 = 0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)                       
                        {
                            double e_case, e_control; //for x2
                            e_case = (t_case[i, j, k, l] + t_control[i, j, k, l] + TO_ZERO) * cnt_case;
                            e_control = (t_case[i, j, k, l] + t_control[i, j, k, l] + TO_ZERO) * cnt_control;
                            x2 += (t_case[i, j, k, l] * cnt - e_case) * (t_case[i, j, k, l] * cnt - e_case) / (e_case + TO_ZERO);
                            x2 += (t_control[i, j, k, l] * cnt - e_control) * (t_control[i, j, k, l] * cnt - e_control) / (e_control + TO_ZERO);
                        }           

            if (1 == opt) //X2
                return x2;	//임상섭
            else
                return (-h_xyt + h_type + h_xy);	//임상섭
        }
        private double distance_5(int a, int b, int c, int d, int e, int opt)
        {
            double h_xy = 0, h_xyt = 0;
            double[, , ,,] table = new double[3, 3, 3, 3,3];
            double[, , ,,] t_case = new double[3, 3, 3, 3,3];
            double[, , ,,] t_control = new double[3, 3, 3, 3,3];
            double[] pa_case = new double[3], pa_control = new double[3],
                pb_case = new double[3], pb_control = new double[3],
                pc_case = new double[3], pc_control = new double[3],
                pd_case = new double[3], pd_control = new double[3],
                pe_case = new double[3], pe_control = new double[3];

            double t_ratio;
            int cnt = 0, cnt_case = 0, cnt_control = 0;
            double x2;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            for(int m=0; m<3; m++)
                            {
                                table[i, j, k, l, m] = 0;
                                t_case[i, j, k, l, m] = 0;
                                t_control[i, j, k, l, m] = 0;
                            }

                pa_case[i] = 0;                pa_control[i] = 0;
                pb_case[i] = 0;                pb_control[i] = 0;
                pc_case[i] = 0;                pc_control[i] = 0;
                pd_case[i] = 0;                pd_control[i] = 0;
                pe_case[i] = 0;                pe_control[i] = 0;
            }

            for (int i = 0; i < num_sample; i++)
            {
                if (snp[a].genotype[i] - '0' >= 0 && snp[b].genotype[i] - '0' >= 0 && snp[c].genotype[i] - '0' >= 0 && snp[d].genotype[i] - '0' >= 0 && snp[e].genotype[i] - '0' >= 0)
                {
                    table[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0', snp[d].genotype[i] - '0', snp[e].genotype[i] - '0']++;
                    cnt++;
                    if (i < num_case)
                    {
                        t_case[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0', snp[d].genotype[i] - '0', snp[e].genotype[i] - '0']++;
                        cnt_case++;
                    }
                    else
                    {
                        t_control[snp[a].genotype[i] - '0', snp[b].genotype[i] - '0', snp[c].genotype[i] - '0', snp[d].genotype[i] - '0', snp[e].genotype[i] - '0']++;
                        cnt_control++;
                    }
                }
            }

            t_ratio = (double)num_case / num_control;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            for(int m=0; m<3; m++)
                            {
                                table[i, j, k, l,m] /= cnt;
                                if (table[i, j, k, l, m] > TO_ZERO)
                                    h_xy += -table[i, j, k, l, m] * Math.Log(table[i, j, k, l, m]) / Math.Log(2.0);

                                t_case[i, j, k, l, m] /= cnt;
                                t_control[i, j, k, l, m] /= cnt;
                                if (t_case[i, j, k, l, m] > TO_ZERO)
                                    h_xyt += -t_case[i, j, k, l, m] * Math.Log(t_case[i, j, k, l, m]) / Math.Log(2.0);
                                if (t_control[i, j, k, l, m] > TO_ZERO)
                                    h_xyt += -t_control[i, j, k, l, m] * Math.Log(t_control[i, j, k, l, m]) / Math.Log(2.0);
                            }


            x2 = 0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            for(int m=0; m<3; m++)
                            {
                                double e_case, e_control; //for x2
                                e_case = (t_case[i, j, k, l, m] + t_control[i, j, k, l, m] + TO_ZERO) * cnt_case;
                                e_control = (t_case[i, j, k, l, m] + t_control[i, j, k, l, m] + TO_ZERO) * cnt_control;
                                x2 += (t_case[i, j, k, l, m] * cnt - e_case) * (t_case[i, j, k, l, m] * cnt - e_case) / (e_case + TO_ZERO);
                                x2 += (t_control[i, j, k, l, m] * cnt - e_control) * (t_control[i, j, k, l, m] * cnt - e_control) / (e_control + TO_ZERO);
                            }

            if (1 == opt) //X2
                return x2;	//임상섭
            else
                return (-h_xyt + h_type + h_xy);	//임상섭
        }

        private double combin(int n, int c)
        {
            double combi = 1;
            for(int i=0; i<c; i++){
                combi *= (n - i);
                combi /= (i + 1);
            }
            return combi;
        }

        private void k_means(int k_num, int opt, int iteration)
        {
	        int cnt, changed, min_index=0;
	        double temp_distance=0, min_distance;
            label9.Text = k_num.ToString() + "-locus clustering";
            label9.Refresh();

            Random rand = new Random();
            Array.Resize<CENTROID>(ref centroid, k_num);
	        //Initialization
	        for(int i=0; i<k_num; i++){
                centroid[i] = new CENTROID();
		        centroid[i].index = rand.Next(0, num_SNP);
		        centroid[i].entropy = snp[i].entropy;
		        centroid[i].ave = 0;
	        }
	        
	        changed = 1;
	        for(int t=0; t<iteration && changed >0 ; t++)
	        {
                progressBar1.Value = (int)(t+1) * 100 / iteration;
                //System.Threading.Thread.Sleep(1);
		        changed = 0;
		        //Finding the closest centroid
		        for(int k=0;k<k_num; k++)
			        centroid[k].ave = 0;
		        for(int j=0; j<num_SNP; j++)
		        {
			        min_index = 0;	min_distance = 10000;
			        for(int k=0; k<k_num; k++)
			        {
                        temp_distance = distance_2(centroid[k].index, j, opt);
                                                    
				        if(min_distance > temp_distance)
				        {
					        min_index = k;
					        min_distance = temp_distance;
				        }
					
			        }
			        if(snp[j].cluster != min_index)
			        {
				        changed++;
				        snp[j].cluster = min_index;				
			        }
			        centroid[min_index].ave += min_distance;
		        }
		
		        //reclaculating centroids		
		        for(int m=0; m<k_num; m++)
		        {
                    double[,] genotype = new double[num_sample, 3];
			        cnt = 0;
			        for(int j=0; j<num_sample ; j++)
				        genotype[j,0] = genotype[j,1] = genotype[j,2] = 0;
			
			        for(int j=0; j<num_SNP; j++)
			        {
				        if(m == snp[j].cluster)
				        {	
					        for(int k=0; k<num_sample; k++)
					        {
						        if(snp[j].genotype[k]-'0' >= 0)
							        genotype[k,snp[j].genotype[k]-'0']++;
					        }
					        cnt++;
				        }
			        }
			        for(int j=0; j<num_sample; j++)
				        for(int k=0; k<3; k++)
				        {
					        genotype[j,k] /= cnt;
					        if(genotype[j,k] < TO_ZERO)
						        genotype[j,k] = TO_ZERO;
				        }
			
			        //대표 sNP 선정
			        cnt = 0; centroid[m].entropy = snp[centroid[m].index].entropy;
			        for(int j=0; j<num_SNP; j++)
			        {
				        if(m == snp[j].cluster)
				        {
                            if(2 == k_num)
                                snp[j].sum = distance_2(centroid[m].index, j, opt);
                            
					        if(0 == cnt || snp[j].sum < snp[min_index].sum)
						        min_index = j;
			
					        cnt++;
				        }
			        }			       
                    
                    centroid[m].index = min_index;
			        centroid[m].ave /= cnt;
			        centroid[m].entropy = snp[min_index].entropy;
			        centroid[m].id = snp[min_index].id;
		        }
	        }
            /*
            for(int i=0; i<k_num; i++)
            {
                double max_sum = 0;
                int max_index = 0;
                for(int j=0; j<num_SNP; j++)
                {
                    if(i == snp[j].cluster){
                        double sum = 0;
                        for(int k = 0; k<k_num; k++)
                            if(2 == k_num)
                            {
                                if(0 == opt)
                                    sum += distance_2(centroid[k].index, j, 0);
                                else if(1 == opt)
                                    sum += distance_2(centroid[k].index, j, 1);
                            }

                        if(max_sum < sum)
                        {
                            max_sum = sum;
                            max_index = j;
                        }
                    }
                }
                centroid[i].id = max_index;
                //printf("(%d-%d)\n", i, centroid[i].id);
             
            }*/
        }

        private void check_set(int k_num, int opt){
            label9.Text = k_num.ToString() + "-locus candidate selection";
            label9.Refresh();
	        //각 SNP의 센트로이드 거리 합 계산하기
	        for(int i=0; i<num_SNP; i++)
	        {
		        snp[i].sum = 0;
		        for(int j=0; j<k_num; j++)
                    if(j != snp[i].cluster)
                        snp[i].sum += distance_2(centroid[j].index, i, opt);
             }

	        //대표 SNP집합 -1로 초기화
	        for(int i=0; i<k_num ; i++)
		        for(int j=0; j<set_size; j++)
			        centroid[i].set[j] = -1;

	        //대표 SNP 집합 선택하기
	        for(int i=0; i<k_num ;i++){
		        for(int j=0; j<num_SNP; j++){
			        if(snp[j].cluster == i){
			             for(int k=0; k<set_size; k++){
					        if(0 > centroid[i].set[k]){	//비어있을 경우
						        centroid[i].set[k] = j;
						        break;
					        }
					        else if(snp[centroid[i].set[k]].sum <= snp[j].sum){
						        for(int l=set_size-1; l>k; l--)
							        centroid[i].set[l] = centroid[i].set[l-1];	//한칸씩 뒤로 밀기
						        centroid[i].set[k] = j;
						        break;
					        }
				        }
			        }
		        }
	        }	        
        }

        private void single_locus()
        {
            label9.Text = "Single-locus";
            label9.Refresh();
            Array.Sort(snp, snpCompare);
            FileStream write_file = new FileStream(textBox1.Text+"_single.txt", FileMode.Create);
            StreamWriter outStream = new StreamWriter(write_file, System.Text.Encoding.Default);
            var d = Gamma.WithShapeScale(1, 1.0/(num_sample*LN2));   //alpha, beta(1/(n*ln2))
            //var d_chi = ChiSquared.CDF(2
            
            string temp;
            temp = "ID\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tmaf\n";
            outStream.WriteLine(temp);
            for(int i=0; i<num_SNP; i++)
            {
                double d_temp;
                temp = snp[i].id.ToString();
                temp += "\t";
                temp += snp[i].mi.ToString();
                temp += "\t";
                if (0 < snp[i].mi)
                {
                    d_temp = ((1.0 - d.CumulativeDistribution(snp[i].mi)) * num_SNP);
                    if (1 > d_temp)
                        temp += d_temp.ToString();
                    else
                        temp += "1.0";
                }
                else
                    temp += "1.0";
                temp += "\t";
                if (0 < snp[i].mi)
                {
                    /*d_temp =  (1 - ChiSquared.CDF(2, (2 * num_sample * LN2 * snp[i].mi))) * num_SNP;
                    if (1 > d_temp)
                        temp += d_temp.ToString();
                    else
                        temp += "1.0";*/
                    temp += "=chisq.dist.rt(";
                    temp += ((2 * num_sample * LN2 * snp[i].mi)).ToString();
                    temp += "," + (2).ToString() + ")";
                    temp += "* combin(" + num_SNP.ToString() + "," + (1).ToString() + ")";
                }
                else
                    temp += "0.0"; 
                temp += "\t";
                temp += snp[i].maf.ToString();
                //temp += "\n";
                outStream.WriteLine(temp);
            }
            outStream.Close();
            write_file.Close();
            label9.Text = "Single-locus end";
        }

        private void two_locus()
        {
            int k_num = 2, opt =0;
            double d_temp;
            label9.Text = "two-locus";
            label9.Refresh();
            k_means(k_num, opt, 10);
            check_set(k_num, opt);
            

            int[] temp_set = new int[set_size * k_num];
            int temp_cnt = 0;            

            for (int i = 0; i < set_size * k_num; i++)
            {
                temp_set[i] = centroid[i / set_size].set[i % set_size];
                if (temp_set[i] >= 0)
                    temp_cnt++;
            }

            int[,] all_combi = new int[(temp_cnt * (temp_cnt - 1)) / 2, k_num];
            double[,] all_value = new double[(temp_cnt * (temp_cnt - 1)) / 2, 4];
            
            label9.Text = "two-locus best_search";
            label9.Refresh();

            int cnt = 0, total = temp_cnt * (temp_cnt-1)/2;
            {
                for (int i = 0; i < set_size * k_num - 1; i++)
                {
                    progressBar1.Value = (int)cnt * 100 / total;
                    for (int j = i + 1; j < set_size * k_num; j++)
                    {                      
                        if (0 <= temp_set[i] && 0 <= temp_set[j])
                        {
                            all_combi[cnt, 0] = temp_set[i];
                            all_combi[cnt, 1] = temp_set[j];
                            all_value[cnt, 0] = distance_2(temp_set[i], temp_set[j], 0);
                            all_value[cnt, 1] = distance_2(temp_set[i], temp_set[j], 1);
                            cnt++; 
                        }
                    }
                }         

            }
            label9.Text = "two-locus sorting";
            label9.Refresh();
            //sorting////////////
            for(int i=0; i<cnt-1; i++)
                for (int j = i + 1; j < cnt; j++)
                {
                    int i_temp;
                    if (all_value[i, opt] < all_value[j, opt])
                    {
                        i_temp = all_combi[i, 0]; all_combi[i, 0] = all_combi[j, 0]; all_combi[j, 0] = i_temp;
                        i_temp = all_combi[i, 1]; all_combi[i, 1] = all_combi[j, 1]; all_combi[j, 1] = i_temp;
                        d_temp = all_value[i, 0]; all_value[i, 0] = all_value[j, 0]; all_value[j, 0] = d_temp;
                        d_temp = all_value[i, 1]; all_value[i, 1] = all_value[j, 1]; all_value[j, 1] = d_temp;
                    }
                }

            var d = Gamma.WithShapeScale((Math.Pow(3, k_num) - 1)/2, 1.0 / (num_sample * LN2));   //alpha, beta(1/(n*ln2)) 
            FileStream write_file = new FileStream(textBox1.Text + "_two_all.txt", FileMode.Create);
            StreamWriter outStream = new StreamWriter(write_file, System.Text.Encoding.Default);
            string s_temp;
            s_temp = "rank\tSNP1\tSNP2\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tboolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream.WriteLine(s_temp);
            int BE_size = Convert.ToInt32(textBox8.Text);
            for (int i = 0; i < cnt; i++)            {
                
                s_temp = (i+1).ToString() +"\t"+ snp[all_combi[i,0]].id.ToString() + "\t";
                s_temp += snp[all_combi[i,1]].id.ToString() + "\t";
                s_temp += all_value[i,0].ToString() + "\t";

                if (0 < all_value[i, 0])
                {
                    d_temp = (1.0 - d.CumulativeDistribution(all_value[i, 0])) * combin(num_SNP, k_num);
                    if (1 > d_temp)
                        s_temp += d_temp.ToString();
                    else
                        s_temp += "1.0";                    
                }
                else
                    s_temp += "1.0";
                s_temp += "\t";
                if (0 < all_value[i, 0])
                {                    
                    s_temp += "=chisq.dist.rt(";
                    s_temp += (2 * num_sample * LN2 * all_value[i, 0]).ToString();
                    s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                    s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                }
                else
                    s_temp += "0.0";

                if (i < BE_size)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


                        for (int k = 0; k < num_sample; k++)
                        {
                            //label12.Text = "for test";
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
                            value = be.run(k_num);
                            //label12.Text = label12.Text + value.ToString();
                            //label12.Refresh();
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                        }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;

                    }
                    s_temp += best.s_print(k_num);
                }
                
                outStream.WriteLine(s_temp);
            }
            outStream.Close();
            write_file.Close();
            
            label9.Text = "two-locus making boolean expression";
            label9.Refresh();
            FileStream write_file1 = new FileStream(textBox1.Text + "_two_deduplicated.txt", FileMode.Create);
            StreamWriter outStream1 = new StreamWriter(write_file1, System.Text.Encoding.Default);
            s_temp = "rank\tSNP1\tSNP2\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_epression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream1.WriteLine(s_temp);

            for (int i = 0; i < cnt; i++)
            {
                progressBar1.Value = (int)(i + 1) * 100 / cnt;
                int d_cnt = 0;
                for(int j=0; j<set_size*k_num;j++)
                    for(int k =0; k<k_num; k++)
                        if(temp_set[j] == all_combi[i, k])
                            d_cnt++;

                label9.Text = "two-locus making boolean expression1";
                label9.Refresh();
                if (d_cnt == k_num)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;
                                                
                        
                        for (int k = 0; k < num_sample; k++)
                        {
                            //label12.Text = "for test";
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
                            value = be.run(k_num);
                            //label12.Text = label12.Text + value.ToString();
                            //label12.Refresh();
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                         }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        //s_temp = be.mi.ToString();
                        //outStream1.WriteLine(s_temp);
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;
                        
                    }

                    s_temp = (i+1).ToString() +"\t" + snp[all_combi[i, 0]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i, 1]].id.ToString() + "\t";
                    s_temp += all_value[i, 0].ToString() + "\t";

                    if (0 < all_value[i, 0])
                    {
                        d_temp = (1.0 - d.CumulativeDistribution(all_value[i, 0])) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";
                    }
                    else
                        s_temp += "1.0";
                    s_temp += "\t";
                    if (0 < all_value[i, 0])
                    {
                        /*d_temp = (1 - ChiSquared.CDF(Math.Pow(3, k_num) - 1, (2 * num_sample * LN2 * all_value[i, 0]))) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";
                         */
                        s_temp += "=chisq.dist.rt(";
                        s_temp += (2 * num_sample * LN2 * all_value[i, 0]).ToString();
                        s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                        s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                    }
                    else
                        s_temp += "0.0";
                    s_temp += best.s_print(k_num);

                    outStream1.WriteLine(s_temp);
                    for (int j = 0; j < set_size * k_num; j++)
                        for (int k = 0; k < k_num; k++)
                            if (temp_set[j] == all_combi[i, k])
                                temp_set[j] = -1;
                }
            }
            outStream1.Close();
            write_file1.Close();
            label9.Text = "two-locus end";
            label9.Refresh();
        }

        private void three_locus()
        {
            int k_num = 3, opt = 0;
            double d_temp;
            label9.Text = "three-locus";
            label9.Refresh();
            k_means(k_num, opt, 10);
            check_set(k_num, opt);


            int[] temp_set = new int[set_size * k_num];
            int temp_cnt = 0;

            for (int i = 0; i < set_size * k_num; i++)
            {
                temp_set[i] = centroid[i / set_size].set[i % set_size];
                if (temp_set[i] >= 0)
                    temp_cnt++;
            }

            int[,] all_combi = new int[(temp_cnt * (temp_cnt - 1) * (temp_cnt - 2)) / 2/3, k_num];
            double[,] all_value = new double[(temp_cnt * (temp_cnt - 1) * (temp_cnt - 2)) / 2 / 3, 4];

            label9.Text = "three-locus best_search";
            label9.Refresh();

            int cnt = 0, total = (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2)) / 2 / 3;
            for (int i = 0; i < set_size * k_num - 2; i++)
            {
                progressBar1.Value = (int)cnt * 100 / total;
                //System.Threading.Thread.Sleep(1);
                for (int i1 = i + 1; i1 < set_size * k_num-1; i1++)
                    for(int i2 = i1+1;i2<set_size*k_num;i2++)
                        if (0 <= temp_set[i] && 0 <= temp_set[i1] && 0 <= temp_set[i2])
                        {
                            all_combi[cnt, 0] = temp_set[i];
                            all_combi[cnt, 1] = temp_set[i1];
                            all_combi[cnt, 2] = temp_set[i2];
                            all_value[cnt, 0] = distance_3(temp_set[i], temp_set[i1], temp_set[i2], 0);
                            all_value[cnt, 1] = distance_3(temp_set[i], temp_set[i1],temp_set[i2], 1);
                            cnt++;
                        }
            }
            label9.Text = "three-locus sorting";
            label9.Refresh();
            //sorting////////////
            for (int i = 0; i < cnt - 1; i++)
                for (int j = i + 1; j < cnt; j++)
                {
                    int i_temp;
                    if (all_value[i, opt] < all_value[j, opt])
                    {
                        for (int k = 0; k < k_num; k++)
                        { i_temp = all_combi[i, k]; all_combi[i, k] = all_combi[j, k]; all_combi[j, k] = i_temp; }
                        for (int k = 0; k < 4; k++)
                        {d_temp = all_value[i, k]; all_value[i, k] = all_value[j, k]; all_value[j, k] = d_temp;}
                     }
                }

            var d = Gamma.WithShapeScale((Math.Pow(3, k_num) - 1) / 2, 1.0 / (num_sample * LN2));   //alpha, beta(1/(n*ln2)) 

            FileStream write_file = new FileStream(textBox1.Text + "_three_all.txt", FileMode.Create);
            StreamWriter outStream = new StreamWriter(write_file, System.Text.Encoding.Default);
            string s_temp;
            s_temp = "rank\tSNP1\tSNP2\tSNP3\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream.WriteLine(s_temp);
            int BE_size = Convert.ToInt32(textBox8.Text);
            for (int i = 0; i < cnt; i++)
            {
                s_temp = (i + 1).ToString() + "\t" + snp[all_combi[i, 0]].id.ToString() + "\t";
                s_temp += snp[all_combi[i, 1]].id.ToString() + "\t";
                s_temp += snp[all_combi[i, 2]].id.ToString() + "\t";
                s_temp += all_value[i, 0].ToString() + "\t";

                if (0 < all_value[i, 0])
                {
                    d_temp = (1.0 - d.CumulativeDistribution(all_value[i, 0])) * combin(num_SNP, k_num);
                    if (1 > d_temp)
                        s_temp += d_temp.ToString();
                    else
                        s_temp += "1.0";
                }
                else
                    s_temp += "1.0";
                s_temp += "\t";
                if (0 < all_value[i, 0])
                {
                    /*
                    d_temp = (1 - ChiSquared.CDF(Math.Pow(3, k_num) - 1, (2 * num_sample * LN2 * all_value[i, 0]))) * combin(num_SNP, k_num);
                    if (1 > d_temp)
                        s_temp += d_temp.ToString();
                    else
                        s_temp += "1.0";*/
                    s_temp += "=chisq.dist.rt(";
                    s_temp += (2 * num_sample * LN2 * all_value[i, 0]).ToString();
                    s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                    s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                }
                else
                    s_temp += "0.0";

                if (i < BE_size)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


                        for (int k = 0; k < num_sample; k++)
                        {
                            //label12.Text = "for test";
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
                            value = be.run(k_num);
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                        }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        //s_temp = be.mi.ToString();
                        //outStream1.WriteLine(s_temp);
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;

                    }
                    s_temp += best.s_print(k_num);
                }

                outStream.WriteLine(s_temp);
            }
            outStream.Close();
            write_file.Close();

            FileStream write_file1 = new FileStream(textBox1.Text + "_three_deduplicated.txt", FileMode.Create);
            StreamWriter outStream1 = new StreamWriter(write_file1, System.Text.Encoding.Default);
            s_temp = "rank\tSNP1\tSNP2\tSNP3\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream1.WriteLine(s_temp);

            for (int i = 0; i < cnt; i++)
            {
                int d_cnt = 0;
                for (int j = 0; j < set_size * k_num; j++)
                    for (int k = 0; k < k_num; k++)
                        if (temp_set[j] == all_combi[i, k])
                            d_cnt++;


                if (d_cnt == k_num)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


                        for (int k = 0; k < num_sample; k++)
                        {
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
                            value = be.run(k_num);
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                        }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;

                    }
                    s_temp = (i + 1).ToString() + "\t" + snp[all_combi[i, 0]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i, 1]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i, 2]].id.ToString() + "\t";
                    s_temp += all_value[i, 0].ToString() + "\t";
                    //s_temp += all_value[i, 1].ToString();// +"\n";
                    if (0 < all_value[i, 0])
                    {
                        d_temp = (1.0 - d.CumulativeDistribution(all_value[i, 0])) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";
                    }
                    else
                        s_temp += "1.0";
                    s_temp += "\t";
                    if (0 < all_value[i, 0])
                    {
                        /*d_temp = (1 - ChiSquared.CDF(Math.Pow(3, k_num) - 1, (2 * num_sample * LN2 * all_value[i, 0]))) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";
                         */
                        s_temp += "=chisq.dist.rt(";
                        s_temp += (2 * num_sample * LN2 * all_value[i, 0]).ToString();
                        s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                        s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                    }
                    else
                        s_temp += "0.0";

                    s_temp += best.s_print(k_num);
                    outStream1.WriteLine(s_temp);
                    for (int j = 0; j < set_size * k_num; j++)
                        for (int k = 0; k < k_num; k++)
                            if (temp_set[j] == all_combi[i, k])
                                temp_set[j] = -1;
                }
            }
            outStream1.Close();
            write_file1.Close();
            label9.Text = "three-locus end";
            label9.Refresh();
        }

        private void four_locus()
        {
            int k_num = 4, opt = 0;
            double d_temp;
            label9.Text = "four-locus";
            label9.Refresh();
            k_means(k_num, opt, 10);
            check_set(k_num, opt);


            int[] temp_set = new int[set_size * k_num];
            int temp_cnt = 0;

            for (int i = 0; i < set_size * k_num; i++)
            {
                temp_set[i] = centroid[i / set_size].set[i % set_size];
                if (temp_set[i] >= 0)
                    temp_cnt++;
            }

            int[,] all_combi = new int[(temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3)) / 2 / 3/4, k_num];
            double[,] all_value = new double[(temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3)) / 2 / 3/4, 4];

            label9.Text = "four-locus best_search";
            label9.Refresh();

            int cnt = 0, total = (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3)) / 2 / 3/4;
            for (int i = 0; i < set_size * k_num - 3; i++)
            {
                progressBar1.Value = (int)cnt * 100 / total;
                //System.Threading.Thread.Sleep(1);
                for (int i1 = i + 1; i1 < set_size * k_num - 2; i1++)
                    for (int i2 = i1 + 1; i2 < set_size * k_num-1; i2++)
                        for (int i3 = i2 + 1; i3 < set_size * k_num; i3++)
                            if (0 <= temp_set[i] && 0 <= temp_set[i1] && 0 <= temp_set[i2] && 0 <= temp_set[i3])
                            {
                                all_combi[cnt, 0] = temp_set[i];
                                all_combi[cnt, 1] = temp_set[i1];
                                all_combi[cnt, 2] = temp_set[i2];
                                all_combi[cnt, 3] = temp_set[i3];
                                all_value[cnt, 0] = distance_4(temp_set[i], temp_set[i1], temp_set[i2], temp_set[i3], 0);
                                all_value[cnt, 1] = distance_4(temp_set[i], temp_set[i1], temp_set[i2], temp_set[i3], 1);
                                cnt++;
                            }
            }
            label9.Text = "four-locus sorting";
            label9.Refresh();
            //sorting////////////
            for (int i = 0; i < cnt - 1; i++)
                for (int j = i + 1; j < cnt; j++)
                {
                    int i_temp;
                    if (all_value[i, opt] < all_value[j, opt])
                    {
                        for (int k = 0; k < k_num; k++)
                        { i_temp = all_combi[i, k]; all_combi[i, k] = all_combi[j, k]; all_combi[j, k] = i_temp; }
                        for (int k = 0; k < 4; k++)
                        {d_temp = all_value[i, k]; all_value[i, k] = all_value[j, k]; all_value[j, k] = d_temp;}
                    }
                }

            var d = Gamma.WithShapeScale((Math.Pow(3, k_num) - 1) / 2, 1.0 / (num_sample * LN2));   //alpha, beta(1/(n*ln2)) 
            FileStream write_file = new FileStream(textBox1.Text + "_four_all.txt", FileMode.Create);
            StreamWriter outStream = new StreamWriter(write_file, System.Text.Encoding.Default);
            string s_temp;
            s_temp = "rank\tSNP1\tSNP2\tSNP3\tSNP4\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream.WriteLine(s_temp);
            int BE_size = Convert.ToInt32(textBox8.Text);
            for (int i = 0; i < cnt; i++)
            {
                s_temp = (i + 1).ToString() + "\t" + snp[all_combi[i, 0]].id.ToString() + "\t";
                s_temp += snp[all_combi[i, 1]].id.ToString() + "\t";
                s_temp += snp[all_combi[i, 2]].id.ToString() + "\t";
                s_temp += snp[all_combi[i, 3]].id.ToString() + "\t";
                s_temp += all_value[i, 0].ToString() + "\t";
                //s_temp += all_value[i, 1].ToString();// +"\n";
                if (0 < all_value[i, 0])
                {
                    d_temp = (1.0 - d.CumulativeDistribution(all_value[i, 0])) * combin(num_SNP, k_num);
                    if (1 > d_temp)
                        s_temp += d_temp.ToString();
                    else
                        s_temp += "1.0";
                }
                else
                    s_temp += "1.0";
                s_temp += "\t";
                if (0 < all_value[i, 0])
                {
                    /*
                    d_temp = (1 - ChiSquared.CDF(Math.Pow(3, k_num) - 1, (2 * num_sample * LN2 * all_value[i, 0]))) * combin(num_SNP, k_num);
                    if (1 > d_temp)
                        s_temp += d_temp.ToString();
                    else
                        s_temp += "1.0";*/
                    s_temp += "=chisq.dist.rt(";
                    s_temp += (2 * num_sample * LN2 * all_value[i, 0]).ToString();
                    s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                    s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                }
                else
                    s_temp += "0.0";

                if (i < BE_size)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


                        for (int k = 0; k < num_sample; k++)
                        {
                            //label12.Text = "for test";
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
                            value = be.run(k_num);
                            //label12.Text = label12.Text + value.ToString();
                            //label12.Refresh();
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                        }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        //s_temp = be.mi.ToString();
                        //outStream1.WriteLine(s_temp);
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;

                    }
                    s_temp += best.s_print(k_num);
                }
                outStream.WriteLine(s_temp);
            }
            outStream.Close();
            write_file.Close();

            FileStream write_file1 = new FileStream(textBox1.Text + "_four_deduplicated.txt", FileMode.Create);
            StreamWriter outStream1 = new StreamWriter(write_file1, System.Text.Encoding.Default);
            s_temp = "rank\tSNP1\tSNP2\tSNP3\tSNP4\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream1.WriteLine(s_temp);

            for (int i = 0; i < cnt; i++)
            {
                int d_cnt = 0;
                for (int j = 0; j < set_size * k_num; j++)
                    for (int k = 0; k < k_num; k++)
                        if (temp_set[j] == all_combi[i, k])
                            d_cnt++;

                if (d_cnt == k_num)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


                        for (int k = 0; k < num_sample; k++)
                        {
                            //label12.Text = "for test";
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
                            value = be.run(k_num);
                            //label12.Text = label12.Text + value.ToString();
                            //label12.Refresh();
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                        }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;

                    }
                    s_temp = (i + 1).ToString() + "\t" + snp[all_combi[i, 0]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i, 1]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i, 2]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i, 3]].id.ToString() + "\t";
                    s_temp += all_value[i, 0].ToString() + "\t";
                    //s_temp += all_value[i, 1].ToString();// +"\n";
                    if (0 < all_value[i, 0])
                    {
                        d_temp = (1.0 - d.CumulativeDistribution(all_value[i, 0])) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";
                    }
                    else
                        s_temp += "1.0";
                    s_temp += "\t";
                    if (0 < all_value[i, 0])
                    {
                        /*d_temp = (1 - ChiSquared.CDF(Math.Pow(3, k_num) - 1, (2 * num_sample * LN2 * all_value[i, 0]))) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";*/
                        s_temp += "=chisq.dist.rt(";
                        s_temp += (2 * num_sample * LN2 * all_value[i, 0]).ToString();
                        s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                        s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                    }
                    else
                        s_temp += "0.0";
                    s_temp += best.s_print(k_num);
                    outStream1.WriteLine(s_temp);
                    for (int j = 0; j < set_size * k_num; j++)
                        for (int k = 0; k < k_num; k++)
                            if (temp_set[j] == all_combi[i, k])
                                temp_set[j] = -1;
                }
            }
            outStream1.Close();
            write_file1.Close();
            label9.Text = "four-locus end";
            label9.Refresh();
        }

        private void five_locus()
        {
            int k_num = 5, opt = 0;
            
            double d_temp;
            label9.Text = "five-locus";
            label9.Refresh();
            k_means(k_num, opt, 10);
            check_set(k_num, opt);


            int[] temp_set = new int[set_size * k_num];
            int temp_cnt = 0;

            for (int i = 0; i < set_size * k_num; i++)
            {
                temp_set[i] = centroid[i / set_size].set[i % set_size];
                if (temp_set[i] >= 0)
                    temp_cnt++;
            }

            COMBINATION[] all_combi = new COMBINATION[(temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3) * (temp_cnt - 4)) / 2 / 3 / 4 / 5];          

            label9.Text = "five-locus best_search";
            label9.Refresh();

            int cnt = 0, total = (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3) * (temp_cnt - 4)) / 2 / 3 / 4/5;
            for (int i = 0; i < set_size * k_num - 4; i++)
            {
                progressBar1.Value = (int)cnt * 100 / total;
                //System.Threading.Thread.Sleep(1);
                for (int i1 = i + 1; i1 < set_size * k_num - 3; i1++)
                    for (int i2 = i1 + 1; i2 < set_size * k_num - 2; i2++)
                        for (int i3 = i2 + 1; i3 < set_size * k_num-1; i3++)
                            for (int i4 = i3 + 1; i4 < set_size * k_num; i4++)
                                if (0 <= temp_set[i] && 0 <= temp_set[i1] && 0 <= temp_set[i2] && 0 <= temp_set[i3] && 0 <= temp_set[i4])
                            {
                                all_combi[cnt] = new COMBINATION();
                                //all_combi[cnt].snp = new int[5];
                                //all_combi[cnt].value = new double[4];
                                all_combi[cnt].snp[0] = temp_set[i];
                                all_combi[cnt].snp[1] = temp_set[i1];
                                all_combi[cnt].snp[2] = temp_set[i2];
                                all_combi[cnt].snp[3] = temp_set[i3];
                                all_combi[cnt].snp[4] = temp_set[i4];
                                all_combi[cnt].value[0] = distance_5(temp_set[i], temp_set[i1], temp_set[i2], temp_set[i3], temp_set[i4], 0);
                                all_combi[cnt].value[1] = distance_5(temp_set[i], temp_set[i1], temp_set[i2], temp_set[i3], temp_set[i4], 1);
                                cnt++;
                            }
            }
            label9.Text = "five-locus sorting";
            label9.Refresh();
            if (0 == opt)
                Array.Sort(all_combi, combiCompare_0);
            else if(1 == opt)
                Array.Sort(all_combi, combiCompare_1);
            else if(2 == opt)
                Array.Sort(all_combi, combiCompare_2);
            else
                Array.Sort(all_combi, combiCompare_3);

            var d = Gamma.WithShapeScale((Math.Pow(3, k_num) - 1) / 2, 1.0 / (num_sample * LN2));   //alpha, beta(1/(n*ln2)) 
            FileStream write_file = new FileStream(textBox1.Text + "_five_all.txt", FileMode.Create);
            StreamWriter outStream = new StreamWriter(write_file, System.Text.Encoding.Default);
            string s_temp;
            s_temp = "rank\tSNP1\tSNP2\tSNP3\tSNP4\tSNP5\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream.WriteLine(s_temp);
            int BE_size = Convert.ToInt32(textBox8.Text);
            for (int i = 0; i < cnt; i++)
            {
                s_temp = (i + 1).ToString() + "\t" + snp[all_combi[i].snp[0]].id.ToString() + "\t";
                s_temp += snp[all_combi[i].snp[1]].id.ToString() + "\t";
                s_temp += snp[all_combi[i].snp[2]].id.ToString() + "\t";
                s_temp += snp[all_combi[i].snp[3]].id.ToString() + "\t";
                s_temp += snp[all_combi[i].snp[4]].id.ToString() + "\t";
                s_temp += all_combi[i].value[0].ToString() + "\t";
                //s_temp += all_combi[i].value[1].ToString();// +"\n";
                if (0 < all_combi[i].value[0])
                {
                    d_temp = (1.0 - d.CumulativeDistribution(all_combi[i].value[0])) * combin(num_SNP, k_num);
                    if (1 > d_temp)
                        s_temp += d_temp.ToString();
                    else
                        s_temp += "1.0";
                }
                else
                    s_temp += "1.0";
                s_temp += "\t";
                if (0 < all_combi[i].value[0])
                {
                    /*d_temp = (1 - ChiSquared.CDF(Math.Pow(3, k_num) - 1, (2 * num_sample * LN2 * all_combi[i].value[0]))) * combin(num_SNP, k_num);
                    if (1 > d_temp)
                        s_temp += d_temp.ToString();
                    else
                        s_temp += "1.0";*/
                    s_temp += "=chisq.dist.rt(";
                    s_temp += (2 * num_sample * LN2 * all_combi[i].value[0]).ToString();
                    s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                    s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                }
                else
                    s_temp += "0.0";

                if (i < BE_size)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


                        for (int k = 0; k < num_sample; k++)
                        {
                            //label12.Text = "for test";
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i].snp[j]].genotype[k] - '0';
                            value = be.run(k_num);
                            //label12.Text = label12.Text + value.ToString();
                            //label12.Refresh();
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                        }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;

                    }                   
                    s_temp += best.s_print(k_num);                    
                }
                outStream.WriteLine(s_temp);
            }
            outStream.Close();
            write_file.Close();

            FileStream write_file1 = new FileStream(textBox1.Text + "_five_deduplicated.txt", FileMode.Create);
            StreamWriter outStream1 = new StreamWriter(write_file1, System.Text.Encoding.Default);
            s_temp = "rank\tSNP1\tSNP2\tSNP3\tSNP4\tSNP5\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN";
            outStream1.WriteLine(s_temp);

            for (int i = 0; i < cnt; i++)
            {
                int d_cnt = 0;
                for (int j = 0; j < set_size * k_num; j++)
                    for (int k = 0; k < k_num; k++)
                        if (temp_set[j] == all_combi[i].snp[k])
                            d_cnt++;

                if (d_cnt == k_num)
                {
                    BOOLEAN_EXPRESSION best = new BOOLEAN_EXPRESSION(k_num);
                    BOOLEAN_EXPRESSION be = new BOOLEAN_EXPRESSION(k_num);
                    best.mi = 0;

                    int t_cnt = 0;
                    while (true)        //test all boolean expression
                    {
                        t_cnt++;
                        int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


                        for (int k = 0; k < num_sample; k++)
                        {
                            for (int j = 0; j < k_num; j++)
                                be.snp[j] = snp[all_combi[i].snp[j]].genotype[k] - '0';
                            value = be.run(k_num);
                            if (1 == value)
                            {
                                if (k < num_case)
                                    case_1++;
                                else
                                    control_1++;
                            }
                            else
                            {
                                if (k < num_case)
                                    case_0++;
                                else
                                    control_0++;
                            }
                            if (case_0 * control_1 > case_1 * control_0)
                            {
                                be.e_gate = 1;
                                be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
                            }
                            else
                            {
                                be.e_gate = 0;
                                be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
                            }
                        }
                        be.mi = 0;
                        be.mi += ((double)case_0 / num_sample) * Math.Log(TO_ZERO + (double)case_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)case_1 / num_sample) * Math.Log(TO_ZERO + (double)case_1 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_0 / num_sample) * Math.Log(TO_ZERO + (double)control_0 / num_sample) / Math.Log(2.0);
                        be.mi += ((double)control_1 / num_sample) * Math.Log(TO_ZERO + (double)control_1 / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + case_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(control_0 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / Math.Log(2.0);

                        be.mi += -((double)(case_0 + control_0) / num_sample) * Math.Log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / Math.Log(2.0);
                        be.mi += -((double)(case_1 + control_1) / num_sample) * Math.Log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / Math.Log(2.0);
                        if (be.mi > best.mi)
                            best.copy(be, k_num);
                        if (1 == be.next(k_num))
                            break;

                    }
                    s_temp = (i + 1).ToString() + "\t" + snp[all_combi[i].snp[0]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i].snp[1]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i].snp[2]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i].snp[3]].id.ToString() + "\t";
                    s_temp += snp[all_combi[i].snp[4]].id.ToString() + "\t";
                    s_temp += all_combi[i].value[0].ToString() + "\t";
                    //s_temp += all_combi[i].value[1].ToString();// +"\n";
                    if (0 < all_combi[i].value[0])
                    {
                        d_temp = (1.0 - d.CumulativeDistribution(all_combi[i].value[0])) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";
                    }
                    else
                        s_temp += "1.0";
                    s_temp += "\t";
                    if (0 < all_combi[i].value[0])
                    {
                        /*d_temp = (1 - ChiSquared.CDF(Math.Pow(3, k_num) - 1, (2 * num_sample * LN2 * all_combi[i].value[0]))) * combin(num_SNP, k_num);
                        if (1 > d_temp)
                            s_temp += d_temp.ToString();
                        else
                            s_temp += "1.0";*/
                        s_temp += "=chisq.dist.rt(";
                        s_temp += (2 * num_sample * LN2 * all_combi[i].value[0]).ToString();
                        s_temp += "," + (Math.Pow(3, k_num) - 1).ToString() + ")";
                        s_temp += "* combin(" + num_SNP.ToString() + "," + k_num.ToString() + ")";
                    }
                    else
                        s_temp += "1.0";

                    s_temp += best.s_print(k_num);
                    outStream1.WriteLine(s_temp);

                    for (int j = 0; j < set_size * k_num; j++)
                        for (int k = 0; k < k_num; k++)
                            if (temp_set[j] == all_combi[i].snp[k])
                                temp_set[j] = -1;
                }
            }
            outStream1.Close();
            write_file1.Close();
            label9.Text = "five-locus end";
            label9.Refresh();
        }

        private void run()
        {
            int cnt = 0;      

            loadData();
            if (true == checkBox1.Checked)
            {
                single_locus();
                cnt++;
            }
            if (true == checkBox2.Checked)
            {
                two_locus();
                cnt++;
            }
            if (true == checkBox3.Checked)
            {
                three_locus();
                cnt++;
            }
            if (true == checkBox4.Checked)
            {
                four_locus();
                cnt++;
            }
            if (true == checkBox5.Checked)
            {
                five_locus();
                cnt++;
            }
            label9.Text = "All " + cnt + " tests are completed.";
            button3.Enabled = true;
        }

        private void button3_Click(object sender, EventArgs e)
        {
            button3.Enabled = false;
            try
            {
                Thread thread = new Thread(new ThreadStart(run));                
                thread.Start();                
            }catch(ThreadStartException err){
                MessageBox.Show(err.Message, err.Source);
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        } 
    }
}
