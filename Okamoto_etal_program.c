/*LINEAR ACCUMULATOR */
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#define L 505 /* cell no. */

/* parameters */
int gg;
int TRIAL;
double dt,time;
int Np,Ni; /* no. of neurons */
int CONNECTIONpp[L][L];
double connection_rate;
double Gpp,Gpi,Gip;

/* Membrane properties of a pyramidal cell */
double V_resting;
double V_theta;
double V_reset;
double G_L;
double refractory;
double Cp;

/* Membrane properties of an inter-neuron */
double Vi_resting;
double Vi_theta;
double Vi_reset;
double G_i_L;
double refractory_i;
double Ci;

/* Membrane potentials and related variables */
double V1[L],V2[L];
double Vi1[L],Vi2[L];
int AP[L];
int APi[L];
double t_last_AP[L];	/* lapse of time after the last AP */
double t_last_APi[L];	/* lapse of time after the last APi */

/* AMPA & GABA */
double E_AMPA,E_GABA;
double tau_AMPA,tau_GABA;
// Sparse connection
double p_AMPA[L],p_GABA[L];
double p_AMPA_jump,p_GABA_jump;

/* Population activity */
double FEFp,FEFi,FREQp,FREQi;

/* Recurrent */
double Ppp,Ppi,Pip;
// Sparse connection
double Ppp_RECURRENT[L];

/* ADP current */
double P_ADP[L],Ca[L],tau_Ca,K_Ca,n_Ca,del_Ca,E_ADP,G_ADP;

/* External current */
double I_base,I_base_i;
double I_dep,I_dep_value,time_dep_on,time_dep_off,I_hyp,I_hyp_value,I_hyp_on,I_hyp_off;

/* Noise */
double G_excitation[L],G_inhibition[L];
double G_excitation_i[L], G_inhibition_i[L];
double F_E,g_E,F_I,g_I;
double F_i_E,g_i_E,F_i_I,g_i_I;


/* Files */
FILE *fp_FEFp,*fp_FEFi;


/* RANDOM */
double RAND()
{
	double b;
	b=(double)rand()/(double)RAND_MAX;
	return(b);
}


/* initialization1 */
initialization1()
{
	char s[80];
	double RAND();
	double pp;
	int n,nn,m,mm,t;

	dt=0.00001;
	TRIAL=1;
	
	/* No. of pyramidal/inhibitory cells */
	Np=100;
	Ni=(int)(0.2*(double)Np);

	/* Membrane propertues of a pyramidal cell */
	V_resting=0.0-70.0;
	V_theta=0.0-52.0;
	V_reset=0.0-62.0;
	Cp=1.0;// scaled to 1.0
	G_L=(1.0/0.02)*Cp;
	refractory=0.002*2.0; 


	/* Membrane propertues of an inter-neuron */
	Vi_resting=0.0-65.0;
	Vi_theta=0.0-52.0;
	Vi_reset=0.0-60.0;
	Ci=1.0;// scaled to 1.0
	G_i_L=(1.0/0.01)*Ci;
	refractory_i=0.001*2.0;


	/* AMPA & GABA */
	E_AMPA=0.0;
	E_GABA=0.0-80.0;
	tau_AMPA=0.005;
	tau_GABA=0.005;
	p_AMPA_jump=0.50;
	p_GABA_jump=0.50;

	/* ADP */
	tau_Ca=0.07;
	K_Ca=1.0;
	n_Ca=4.0;
	del_Ca=1.3;
	E_ADP=0.0-35.0;
	G_ADP=12.0*Cp;

	/* Synaptic connection (conductance) */
	Gpp=20.0/(double)Np*Cp;	/* p->p */
	Gip=4.0/(double)Ni*Ci;		/* p->i */
	Gpi=20.0/(double)Np*Cp;	/* i->p */
	connection_rate=0.1;
		
	/* Noise in a pyramidal cell */
	F_E=(dt/1.0)*(1800.0);
	F_I=F_E;
	g_E=1.4*Cp;
	g_I=g_E*Cp;

	/* Noise in an interneuron */
	F_i_E=(dt/1.0)*(1800.0);
	F_i_I=F_i_E;
	g_i_E=1.4*Ci;
	g_i_I=g_i_E*Ci;

	/* Depolarization/hyperpolarization during the delay */
	// Depolarization
	I_dep_value=50.0*Cp;
	time_dep_on=0.0;
	time_dep_off=100.0;
	// Hyperpolarization
	I_hyp_value=(0.0-0.0*200.0)*Ci;
	I_hyp_on=4.0;
	I_hyp_off=4.2;

	/* Base current */
	I_base=400.0*Cp;
	I_base_i=1000.0*Ci;

	/* Connection matrix */
	for(n=1;n<=Np;++n)
	{
		for(nn=1;nn<=Np;++nn)
		{
			CONNECTIONpp[n][nn]=0;
			pp=RAND();
			if(pp<connection_rate)
			{
				CONNECTIONpp[n][nn]=1;
			}
		}
	}
	for(n=1;n<=Np;++n)
	{
		CONNECTIONpp[n][n]=0;
	}

}


/* initialization2 */
initialization2()
{
	char s[80];
	double RAND();
	int n;
	double pp;

	for(n=1;n<=Np;++n)
	{
		V2[n]=V_resting;
		t_last_AP[n]=10000.0;
		p_AMPA[n]=0.0;
		Ca[n]=0.0;
		G_excitation[n]=0.0;
		G_inhibition[n]=0.0;
	}

	for(n=1;n<=Ni;++n)
	{
		Vi2[n]=Vi_resting;
		t_last_APi[n]=10000.0;
		p_GABA[n]=0.0;
		G_excitation_i[n]=0.0;
		G_inhibition_i[n]=0.0;
	}

	FREQp=0.0;
	FREQi=0.0;

}

rewrite()
{
	int n,nn,m,mm;
	double RAND();

	for(n=1;n<=Np;++n)
	{
		V1[n]=V2[n];
		t_last_AP[n]=t_last_AP[n]+dt;
	}

	for(n=1;n<=Ni;++n)
	{
		Vi1[n]=Vi2[n];
		t_last_APi[n]=t_last_APi[n]+dt;
	}

	Ppp=0.0;

	for(n=1;n<=Np;++n)
	{
/*
		Ppp_RECURRENT[n]=0.0;
		for(nn=1;nn<=Np;++nn)
		{
			Ppp_RECURRENT[n]=Ppp_RECURRENT[n]
			+(double)CONNECTIONpp[n][nn]*p_AMPA[nn]*(1.0/connection_rate);
		}
*/
		Ppp=Ppp+p_AMPA[n];
		P_ADP[n]=pow(Ca[n],n_Ca)/(pow(Ca[n],n_Ca)+pow(K_Ca,n_Ca));
			
	}

	Pip=0.0;
	for(n=1;n<=Np;++n)
	{
		Pip=Pip+p_AMPA[n];
	}

	Ppi=0.0;
	for(n=1;n<=Ni;++n)
	{
		Ppi=Ppi+p_GABA[n];
	}

}


/* update */
update()
{
	int n,m,nn,mm;
	char s[80];
	double gauss();
	double RAND();
	double pp,pp_e,pp_i;

	/* Population activity */
	FREQp=FREQp*(1.0-dt/tau_AMPA);
	FREQi=FREQi*(1.0-dt/tau_GABA);
	FEFp=FREQp/(double)Np;
	FEFi=FREQi/(double)Ni;


	/* Excitatory cells */
	for(n=1;n<=Np;++n)
	{
		p_AMPA[n]=p_AMPA[n]*(1.0-dt/tau_AMPA);
		Ca[n]=Ca[n]*(1.0-dt/tau_Ca);
	
		pp_e=RAND();
		if(pp_e<F_E)
		{
			G_excitation[n]=G_excitation[n]+g_E;
		}
		G_excitation[n]=G_excitation[n]*(1.0-dt/tau_AMPA);
		pp_i=RAND();
		if(pp_i<F_I)
		{
			G_inhibition[n]=G_inhibition[n]+g_I;
		}
		G_inhibition[n]=G_inhibition[n]*(1.0-dt/tau_GABA);

		if(t_last_AP[n]>refractory)
		{
			V2[n]=V1[n]+dt*
			(
			(-1.0)*G_L*(V1[n]-V_resting)/Cp
			/*
			+Gpp*Ppp_RECURRENT[n]*(0.0-V1[n])/Cp
			*/
			+Gpp*(Ppp-p_AMPA[n])*(E_AMPA-V1[n])/Cp
			+Gpi*Ppi*(E_GABA-V1[n])/Cp
			+G_ADP*(E_ADP-V1[n])*P_ADP[n]/Cp
			+I_base/Cp
			+G_excitation[n]*(E_AMPA-V1[n])/Cp
			+G_inhibition[n]*(E_GABA-V1[n])/Cp
			+I_dep/Cp
			+I_hyp/Cp
			);
		}
		else
		{
			V2[n]=V1[n];
		}

		AP[n]=0;
		if(V2[n]>V_theta)
		{
			AP[n]=1;
			p_AMPA[n]=p_AMPA[n]+p_AMPA_jump*(1.0-p_AMPA[n]);
			t_last_AP[n]=0.0;
			/**/	
			FREQp=FREQp+1.0;
			V2[n]=V_reset;
			Ca[n]=Ca[n]+del_Ca;
		}
	}


	/* Inhibitory cells */	
	for(n=1;n<=Ni;++n)
	{
		p_GABA[n]=p_GABA[n]*(1.0-dt/tau_GABA);
		pp_e=RAND();
		if(pp_e<F_i_E)
		{
			G_excitation_i[n]=G_excitation_i[n]+g_i_E;
		}
		G_excitation_i[n]=G_excitation_i[n]*(1.0-dt/tau_AMPA);
		pp_i=RAND();
		if(pp_i<F_i_I)
		{
			G_inhibition_i[n]=G_inhibition_i[n]+g_i_I;
		}
		G_inhibition_i[n]=G_inhibition_i[n]*(1.0-dt/tau_GABA);

		if(t_last_APi[n]>refractory_i)
		{
			Vi2[n]=Vi1[n]+dt*
			(
			(-1.0)*G_i_L*(Vi1[n]-Vi_resting)/Ci
			+Gip*Pip*(E_AMPA-Vi1[n])/Ci
			+I_base_i/Ci
			+G_excitation_i[n]*(E_AMPA-Vi1[n])/Ci
			+G_inhibition_i[n]*(E_GABA-Vi1[n])/Ci
			);
		}
		else
		{
			Vi2[n]=Vi1[n];
		}

		APi[n]=0;
		if(Vi2[n]>Vi_theta)
		{
			APi[n]=1;
			p_GABA[n]=p_GABA[n]+p_GABA_jump*(1.0-p_GABA[n]);
			FREQi=FREQi+1.0;
			Vi2[n]=Vi_reset;
			t_last_APi[n]=0.0;
		}
	}

}

/* insert blanck line */
blanck_line()
{
	fprintf(fp_FEFp,"\n");
	fprintf(fp_FEFi,"\n");
}


/* file_print1 */
file_fprint1()
{
	fprintf(fp_FEFp,"%lf %lf\n",time,FEFp);
	fprintf(fp_FEFi,"%lf %lf\n",time,FEFi);
}


I_on_off()
{
	int n,nn,m,mm;
	double RAND();

		I_dep=0.0;
		if((time>time_dep_on)&&(time<time_dep_off))
		{
			I_dep=I_dep_value;
		}

		I_hyp=0.0;
		if((time>I_hyp_on)&&(time<I_hyp_off))
		{
			I_hyp=I_hyp_value;
		}
}



/* main program */ main() 
{
	char s[80];
	int R_SEED;

	/* file open */
	if((fp_FEFp=fopen("f_FEFp.out","a"))==NULL)
	{printf("File f_FEFp.out not exist\n");}
	if((fp_FEFi=fopen("f_FEFi.out","a"))==NULL)
	{printf("File f_FEFi.out not exist\n");}

	/* set a seed for random variables */
	printf("R_SEED:");
	sscanf(gets(s),"%d",&R_SEED);
	printf("\n");
	srand(R_SEED);

	initialization1();

	for(gg=1;gg<=TRIAL;++gg)
	{
		initialization2();
		time=0.0-2.0;

		blanck_line();
	//	system("del f_scatter.out");

		while(time<3.0)
		{
			time=time+dt;
			rewrite();
			I_on_off();
			update();
			if(((int)(time/dt)%50)==0)
			{
				file_fprint1();
			}
		}


		if((gg%1)==0)
		{
			printf("trial=%d\n",gg);
		}


	}

	/* file close */
	fclose(fp_FEFp);
	fclose(fp_FEFi);

	printf("Program terminated\n");
}

