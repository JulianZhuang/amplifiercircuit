//File Name: Amplifier_circuit_F2020_E50
//Author: Qing Zhuang (Chun Wei Zhuang) 
//Description: Calculate the output volage through numeric and analystic method.
//Last Changed: December 5, 2020
//Running Result:
/*
Please enter file name (do not include the extention): output
Please input the size of dt:  1
Please input the size of dt:  0.5
Please input the size of dt:  0.1
Please input the size of dt:  0.05
Please input the size of dt:  0.01

 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 <<                END OF DISLIN / VERSION 11.3.1                <<
 <<  Date    : 05.12.2020  Time    : 20:45:00  Pageformat: DA4L  <<
 <<  Vectors : 7774        Warnings: 0         Fileformat:  PDF  <<
 <<  Metafile: dislin_6.pdf                                      <<
 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


--------------------------------
Process exited with return value 0
Press any key to continue . . .
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <windows.h>
#include "c:\dislin\dislin.h"
using namespace std;

class Amplifier
{
public:
	Amplifier();

	void input_dt();
	void set_max_analy();
	void set_min_analy();
	void set_max_num();
	void set_min_num();
	int get_n();
		
	float* get_t_array();
	float* get_v_array_analy();
	float* get_v_array_num();

	//v_max
	float get_v_max_analy();
	float get_v_min_analy();
	float get_v_max_num();
	float get_v_min_num();

	//t_max
	float get_t_max_analy();
	float get_t_min_analy();
	float get_t_max_num();
	float get_t_min_num();
	
	
	float v_analytical(float t);
	float v_numeric(float t, float v_0, float v_1); //Using differential equation

	void assign_t_array(float t[]);
	void assign_v_aray_analytical(float t[], float v[]);
	void assign_v_aray_numeric(float t[], float v[]);
	
	void assign_arrays();
	void delete_arrays();

private:
	float v0_num=0., v1_num=0.;
	float dt=0.01;
	int n= 15/dt;
	
	//For return
	float  *t_array, *v_array_analy, *v_array_num;
	float v_max_analy, v_min_analy, v_max_num, v_min_num;
	float t_max_analy, t_min_analy, t_max_num, t_min_num;
};

string name_filename();

/*========================================================================*/
int main()
{
	int ic; //For dislin plot only
	
	ofstream write;
	write.open(name_filename().c_str());
	write << fixed;
	
	int test_num = 1, write_width=15;
	Amplifier *test = new Amplifier[test_num];


	for(int i=0; i<test_num; i++)
	{
		int n = test[i].get_n();
		float *t = test[i].get_t_array();
		float *v_analy = test[i].get_v_array_analy();
		float *v_num = test[i].get_v_array_num();
		
		write << "\n\n\n===================================================================================\n";
		write <<setw(write_width) << "TEST " << i << "\n";
		
		write << "================Maximum & Minimum==================================================\n\n";
		write << setw(write_width) << "t_max_analy: " << test[i].get_t_max_analy() << setw(write_width) 
			<< "v_max_analy: " << test[i].get_v_max_analy() << "\n";
		write << setw(write_width) << "t_min_analy: " << test[i].get_t_min_analy() << setw(write_width) 
			<< "v_min_analy: " << test[i].get_v_min_analy() << "\n";
		write << setw(write_width) << "t_max_num: " << test[i].get_t_max_num() << setw(write_width) 
			<< "v_max_num: " << test[i].get_v_max_num() << "\n";	
		write << setw(write_width) << "t_min_num: " << test[i].get_t_min_num() << setw(write_width) 
			<< "v_min_num: " << test[i].get_v_min_num() << "\n";	
		
		write << "=================t - v table=======================================================\n\n";		
			
		write << setw(write_width) << "t" << setw(write_width) << "v_analy" << setw(write_width) 
			<< "v_num" << "\n\n";
		for(int j=0; j<n; j++)
		{
			if(j%10==0)
			{
			write << setw(write_width) << t[j] << setw(write_width) << v_analy[j] << setw(write_width)
				<< v_num[j] << "\n";
			}	
		}
			
		
	}
	
	//Output Summary
	write << "\n\n=================Summary===========================================================\n";
	write << "List of classes and their structure: We only have one class: Amplifier\n\n";
	write << "Two user designed functions include:\n1.Amplifier::input_dt()\n2.name_filename()\n\n";
	write << "Discuss the differences between the two solutions, the effect of delta t, the time step, on the accuracy of the result:\n";
	write << "According to the graph and table, we find that the analytical approach get accuracy faster than the numric one at the same time step.\n";	
	write << "As the dt (time step) get smaller, the accuracy of numeric approach improves and gets closer to analytical one. Similarly for analytical one," 
		<<"\nhowever, it's more about the design of curve function in which how many connected points there are, whose concept is not quite the same as differential equation.\n";
	//***************Plot Configuration******************    
	metafl("pdf");//Creates screen output. To create PDF output use "pdf"  
	disini();  
	pagera();  
	complx();  
	axspos(450,1800);  
	axslen(2200,1200);  
	name("Time (sec)","x");  
	name("Voltages","y");  
	labdig(-1,"x");  
	ticks(9,"x");  
	ticks(10,"y");  
	titlin("Voltage over Time",1);  
	titlin("Analytical vs Numeric",3);  
	ic=intrgb(0.95,0.95,0.95);  
	axsbgd(ic);  
	graf(0.,15.,0.,1.,-1.,20.,0.,1.);  
	setrgb(0.7,0.7,0.7);  
	grid(1,1);  
	color("fore");  
	height(50);  
	title();
	//********************Create Plot********************  
	
	//Analytical
	color("red");
	for(int i=0; i<test_num; i++)
	{
		curve(test[i].get_t_array(),test[i].get_v_array_analy(),test[i].get_n()); //plots
	}
 
	//Numeric
	color("green"); 
	for(int i=0; i<test_num; i++)
	{
		curve(test[i].get_t_array(),test[i].get_v_array_num(),test[i].get_n()); //plots
	}

	disfin(); //Ends plot session. Must include.
	//***************************************************  

	
	//Memory Management
	write.close();
	
	for(int i=0; i<test_num; i++)
	//Delete all inside-class arrays
	{
		test[i].delete_arrays();
	}
	delete [] test;

	return 0;
}
/*========================================================================*/
Amplifier::Amplifier()
{
	input_dt();
	
	t_array = new float[n];
	v_array_analy = new float[n];
	v_array_num = new float[n];	
	
	assign_arrays();

	set_max_analy();
	set_min_analy();
	set_max_num();
	set_min_num();		
}


void Amplifier::input_dt()
{
	float num;
	cout << "Please input the size of dt:  ";
	cin >> num;
	dt = num;
	n = 15/dt;
}


void Amplifier::set_max_analy()
{
	float left, right;
	t_max_analy = t_array[0], v_max_analy = v_array_analy[0];
	
	for(int i=0; i<n-1; i++)
	{
		left = v_array_analy[i];
		right = v_array_analy[i+1];

		if(left<right && v_max_analy<right)
		{
			t_max_analy = t_array[i];
			v_max_analy = right;
		}		
	}
}

void Amplifier::set_min_analy()
{
	float left, right; 
	t_min_analy = t_array[0], v_min_analy = v_array_analy[0];
	
	for(int i=0; i<n-1; i++)
	{
		left = v_array_analy[i];
		right = v_array_analy[i+1];

		if(left>right && v_min_analy>right)
		{
			t_min_analy = t_array[i];
			v_min_analy = right;
		}		
	}
}

void Amplifier::set_max_num()
{
	float left, right; 
	t_max_num = t_array[0], v_max_num = v_array_num[0];
	
	for(int i=0; i<n-1; i++)
	{
		left = v_array_num[i];
		right = v_array_num[i+1];

		if(left<right && v_max_num<right)
		{
			t_max_num = t_array[i];
			v_max_num = right;
		}		
	}
}

void Amplifier::set_min_num()
{
	float left, right;
	t_min_num = t_array[0], v_min_num = v_array_num[0];
	
	for(int i=0; i<n-1; i++)
	{
		left = v_array_num[i];
		right = v_array_num[i+1];

		if(left>right && v_min_num>right)
		{
			t_min_num = t_array[i];
			v_min_num = right;
		}		
	}
}

int Amplifier::get_n() {return n;}	

	
float* Amplifier::get_t_array()
{
	return t_array;
}

float* Amplifier::get_v_array_analy()
{
	return v_array_analy;
}

float* Amplifier::get_v_array_num()
{
	return v_array_num;
}	

//v_max
float Amplifier::get_v_max_analy()
{
	return v_max_analy;
}

float Amplifier::get_v_min_analy()
{
	return v_min_analy;
}

float Amplifier::get_v_max_num()
{
	return v_max_num;
}	

float Amplifier::get_v_min_num()
{
	return v_min_num;
}
//t_max
float Amplifier::get_t_max_analy()
{
	return v_max_analy;
}

float Amplifier::get_t_min_analy()
{
	return t_min_analy;
}

float Amplifier::get_t_max_num()
{
	return t_max_num;
}	

float Amplifier::get_t_min_num()
{
	return t_min_num;
}


float Amplifier::v_analytical(float t)
{
	float v;
	v = 10 - exp(-t)*(10*cos(2*t)+5*sin(2*t));
	return v;
}

float Amplifier::v_numeric(float t, float v_0, float v_1)
//Using differential equation
{
	float v_2;
	v_2 = pow( (pow(1/dt,2) + (2/dt)), -1)*
		(50 - (v_0/pow(dt, 2)) + (v_1*(-5 + (2/pow(dt, 2)+(2/dt)))));
		
	return v_2;
}

void Amplifier::assign_t_array(float t[])
{
	float increment = 0.;
	
	for(int i=0; i<n; i++)
	{
		t[i]=increment;
		increment+=dt;
	}
}

void Amplifier::assign_v_aray_analytical(float t[], float v[])
{	
	for(int i=0; i<n; i++)
	{
		v[i] = v_analytical(t[i]);
	}
}

void Amplifier::assign_v_aray_numeric(float t[], float v[])
{
	float v_old = v0_num, v_new = v1_num;
	v[0] = v0_num; v[1] = v1_num;
	
	for(int i=2; i<n; i++)
	{
		v[i] = v_numeric(t[i], v_old, v_new);
		v_old = v_new;
		v_new = v[i];
	}
}

void Amplifier::assign_arrays()
{
	assign_t_array(t_array);
	assign_v_aray_analytical(t_array, v_array_analy);
	assign_v_aray_numeric(t_array, v_array_num);	
}

void Amplifier::delete_arrays()
{
	delete [] t_array; 
	delete [] v_array_analy; 
	delete [] v_array_num;
}

string name_filename()
{
	string filename;
	cout << "Please enter file name (do not include the extention): ";
	cin >> filename;
	
	filename = filename + ".txt";
	
	return filename;
}
