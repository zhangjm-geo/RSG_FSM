/*------------------------------------------------------------------------
 *  A rotated staggered grid fast sweeping method for tralveltime calculation
 *
 *  Jing Wang and Jianming Zhang, 9.3.2026
 * 
 *  Contanct: zhangjianming@ouc.edu.cn
 *  ----------------------------------------------------------------------*/
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>
#include<time.h>
using namespace std;
float pi=3.1415926;

float min(float a,float b)
{
	if (a>b) return b;
	else return a;
}

int main(int argc, char **argv)
{
	int nx, nz;
	float dx, dz, velocity;
	float s0 = 0.0;
	float G0 = 0.0;
	int x0, z0;
	int method, analy;
	char name[256];
	clock_t start, end;

	char *file_para;
	file_para = argv[1];

	FILE *fp;
	fp = fopen(file_para, "r");
	if(fp == NULL)
	{
		printf("\n==================================================================\n");
		printf(" Cannot open input file %s \n",file_para);
		printf("\n==================================================================\n\n");
		exit(0);
	}
/*============================================read parameters============================================================*/
	char velmodel[256];
	char line[256];

	fgets(line, sizeof(line), fp);
	sscanf(line, "--x grid number: nx %d", &nx);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--z grid number: nz %d", &nz);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--x grid interval (m): dx %f", &dx);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--z grid interval (m): dz %f", &dz);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--x source grid position: x0 %d", &x0);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--z source grid position: z0 %d", &z0);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--velocity model file name: velmodel %s", velmodel);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--method(1_FSM;2_2ndFSM;3_SGFSM): method %d", &method);
	
	fgets(line, sizeof(line), fp);
	sscanf(line, "--analy(1_homo;2_gradient_slow;3_no): analy %d", &analy);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--velocity model: velocity %f", &velocity);

	fgets(line, sizeof(line), fp);
	sscanf(line, "--constant gradient: gradient %f", &G0);
	fclose(fp);

	int i,j;
	float **s=new float*[nx], **temp = new float*[nx], **T_old=new float*[nx], **T_ana=new float*[nx], **T_new=new float*[nx],**T=new float*[nx],**delta=new float*[nx],**eta=new float*[nx],**v=new float*[nx],**vt=new float*[nx];
	for(int i=0;i<nx;i++)
	{
		s[i]=new float[nz]; temp[i] = new float[nz];T_old[i]=new float[nz]; T_ana[i]=new float[nz]; T_new[i]=new float[nz];T[i]=new float[nz];delta[i]=new float[nz];eta[i]=new float[nz];v[i]=new float[nz];vt[i]=new float[nz];
		for(int j=0;j<nz;j++)
		{	s[i][j]=0; temp[i][j] = 0; T_old[i][j]=0; T_ana[i][j]=0; T_new[i][j]=60;T[i][j]=0;delta[i][j]=0;eta[i][j]=0;v[i][j]=0;vt[i][j]=0;}
	}

/*==================================================read velocity model======================================================*/
	if(analy == 1) //creat constant velocity model
	{
		fp = fopen(velmodel, "wb+");
		for(i=0;i<nx;i++)
			for(j=0;j<nz;j++)
			{
				v[i][j]=velocity;
				fwrite(&v[i][j],sizeof(float),1,fp);
			}
		fclose(fp);	
	}
	else if(analy == 2) //creat constant gradient velocity model
	{
		fp = fopen(velmodel, "wb+");
		for(i=0;i<nx;i++)
			for(j=0;j<nz;j++)
			{
				v[i][j]=velocity+sqrt(pow((i-x0)*dx,2)+pow((j-z0)*dz,2))*G0;
				//v[i][j]=velocity+(j-z0)*dz*G0;
				fwrite(&v[i][j],sizeof(float),1,fp);
			}
		fclose(fp);	
	}

	fp = fopen(velmodel, "r");
	for(i=0;i<nx;i++)
		for(j=0;j<nz;j++)
			fread(&v[i][j], sizeof(float), 1, fp);
	fclose(fp);


	for(i=0;i<nx;i++)
	{
		for(j=0;j<nz;j++)   
		{
			s[i][j]=(1.0/v[i][j]);
			if(i == x0 && j == z0) s0 = s[x0][z0];
			T_old[i][j]=22.0;
			T[i][j]=22.0;
		}  
	}   


	T_old[x0][z0]=0.0;

	float T_xmin=33.0; float T_zmin=33.0;
	float tmp;
	int loop=0,Max_loop=1000;
	int sx,sz;
	int xa,xb,za,zb;
	int xa1,xb1,za1,zb1;
	float lnorm = 60.0, tnorm = 1e-20;

	start = clock();
/*===============================================calculate traveltime=========================================================*/
	while(lnorm >= tnorm && loop < Max_loop)
	{
		//i++ j++
		for(i=0;i<nx;i++)
			for(j=0;j<nz;j++)
			{
				xa=i+1;xb=i-1;za=j+1;zb=j-1;
				if(xa==nx) xa=i;
				if(xb==-1) xb=i;
				if(za==nz) za=j;
				if(zb==-1) zb=j;

				//1stFD
				if(method == 1 || method == 2)
				{
					T_xmin=min(T_old[xa][j],T_old[xb][j]);
					T_zmin=min(T_old[i][za],T_old[i][zb]);
					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+s[i][j]*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);

				}


				//RSGFSM
				if(method == 3)
				{
					if(T_old[i][j]>=T_old[xb][zb]&&T_old[i][j]>=T_old[i][zb]&&T_old[i][j]>=T_old[xb][j])
					{
						T_old[i][zb]=min(T_old[xb][zb]+s[xb][zb]*dx,T_old[i][zb]);
						T_old[xb][j]=min(T_old[xb][zb]+s[xb][zb]*dz,T_old[xb][j]);
						tmp=pow(s[xb][zb]*sqrt(dx*dx+dz*dz),2)-pow(T_old[i][zb]-T_old[xb][j],2);
						if(tmp>=0) T_old[i][j]=min(T_old[xb][zb]+sqrt(tmp),T_old[i][j]);
					}
				}


			}

		//i-- j++
		for(i=nx-1;i>=0;i--)
			for(j=0;j<nz;j++)
			{
				xa=i+1;xb=i-1;za=j+1;zb=j-1;
				if(xa==nx) xa=i;
				if(xb==-1) xb=i;
				if(za==nz) za=j;
				if(zb==-1) zb=j;
				
				//1stFD
				if(method == 1 || method == 2)
				{
					T_xmin=min(T_old[xa][j],T_old[xb][j]);
					T_zmin=min(T_old[i][za],T_old[i][zb]);
					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+s[i][j]*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);
					
				}

				//RSGFSM
				if(method == 3)
				{
					if(T_old[i][j]>=T_old[xa][zb]&&T_old[i][j]>=T_old[xa][j]&&T_old[i][j]>=T_old[i][zb])
					{
						T_old[xa][j]=min(T_old[xa][zb]+s[xa][zb]*dx,T_old[xa][j]);
						T_old[i][zb]=min(T_old[xa][zb]+s[xa][zb]*dz,T_old[i][zb]);
						tmp=pow(s[xa][zb]*sqrt(dx*dx+dz*dz),2)-pow(T_old[i][zb]-T_old[xa][j],2);
						if(tmp>=0)
						T_old[i][j]=min(T_old[xa][zb]+sqrt(tmp),T_old[i][j]);	
					}						
				}

			}

		//i-- j--
		for(i=nx-1;i>=0;i--)
			for(j=nz-1;j>=0;j--)
			{
				xa=i+1;xb=i-1;za=j+1;zb=j-1;
				if(xa==nx) xa=i;
				if(xb==-1) xb=i;
				if(za==nz) za=j;
				if(zb==-1) zb=j;

				//1stFD
				if(method == 1 || method == 2)
				{
					T_xmin=min(T_old[xa][j],T_old[xb][j]);
					T_zmin=min(T_old[i][za],T_old[i][zb]);
					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+s[i][j]*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);
				}

				//RSGFSM
				if(method == 3)
				{
					if(T_old[i][j]>=T_old[xa][za]&&T_old[i][j]>=T_old[xa][j]&&T_old[i][j]>=T_old[i][za])
					{
						T_old[xa][j]=min(T_old[xa][za]+s[xa][za]*dx,T_old[xa][j]);
						T_old[i][za]=min(T_old[xa][za]+s[xa][za]*dz,T_old[i][za]);
						tmp=pow(s[xa][za]*sqrt(dx*dx+dz*dz),2)-pow(T_old[i][za]-T_old[xa][j],2);
						if(tmp>=0) T_old[i][j]=min(T_old[xa][za]+sqrt(tmp),T_old[i][j]);
					}					
				}

			}

		//i++ j--
		for(i=0;i<nx;i++)
			for(j=nz-1;j>=0;j--)
			{
				xa=i+1;xb=i-1;za=j+1;zb=j-1;
				if(xa==nx) xa=i;
				if(xb==-1) xb=i;
				if(za==nz) za=j;
				if(zb==-1) zb=j;

				//1stFD
				if(method == 1 || method == 2)
				{
					T_xmin=min(T_old[xa][j],T_old[xb][j]);
					T_zmin=min(T_old[i][za],T_old[i][zb]);
					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+(s[i][j])*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);
					
				}
				
				//RSGFSM
				if(method == 3)
				{
					if(T_old[i][j]>=T_old[xb][za]&&T_old[i][j]>=T_old[xb][j]&&T_old[i][j]>=T_old[i][za])
					{
						T_old[xb][j]=min(T_old[xb][za]+s[xb][za]*dx,T_old[xb][j]);
						T_old[i][za]=min(T_old[xb][za]+s[xb][za]*dz,T_old[i][za]);
						tmp=pow(s[xb][za]*sqrt(dx*dx+dz*dz),2)-pow(T_old[i][za]-T_old[xb][j],2);
						if(tmp>=0) T_old[i][j]=min(T_old[xb][za]+sqrt(tmp),T_old[i][j]);
					}											
				}

			}
		
			lnorm = 0.0;
			for (i=0;i<nx;i++)
			{
				for (j=0;j<nz;j++)
				{							
					lnorm += fabs(T_new[i][j]-T_old[i][j]);							
				}
			}

			for (i=0;i<nx;i++)
			{
				for (j=0;j<nz;j++)
				{							
					T_new[i][j]= T_old[i][j];							
				}
			}									
			loop++;
	}//end of loop

	//2ndFSM
	if(method == 2)
	{   
		float Txf;float Txz;float Tzf; float Tzz;
		float wxf,wxz,wzf,wzz,rxf,rxz,rzf,rzz; float epsi=1e-6;
		loop=0, lnorm = 60;
		while(lnorm >= tnorm && loop < Max_loop)		
		{			
			for(i=0;i<nx;i++)
				for(j=0;j<nz;j++)
				{
					xa=i+1;xb=i-1;za=j+1;zb=j-1;
					if(xa==nx) xa=i;
					if(xb==-1) xb=i;
					if(za==nz) za=j;
					if(zb==-1) zb=j;
					
					xa1=i+2;xb1=i-2;za1=j+2;zb1=j-2;
					if(xa1==nx) xa1=xa;
					if(xb1==-1) xb1=xb;
					if(za1==nz) za1=za;
					if(zb1==-1) zb1=zb;
					if(xa1==nx+1) xa1=i;
					if(xb1==-2) xb1=i;
					if(za1==nz+1) za1=j;
					if(zb1==-2) zb1=j;
				
					rxf=(epsi+pow(T_old[i][j]-2*T_old[xb][j]+T_old[xb1][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rxz=(epsi+pow(T_old[xa1][j]-2*T_old[xa][j]+T_old[i][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rzf=(epsi+pow(T_old[i][j]-2*T_old[i][zb]+T_old[i][zb1],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					rzz=(epsi+pow(T_old[i][za1]-2*T_old[i][za]+T_old[i][j],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					wxf=1.0/(1.0+2.0*pow(rxf,2));
					wxz=1.0/(1.0+2.0*pow(rxz,2));
					wzf=1.0/(1.0+2.0*pow(rzf,2));
					wzz=1.0/(1.0+2.0*pow(rzz,2));
					Txf=(1.0-wxf)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxf*(3.0*T_old[i][j]-4.0*T_old[xb][j]+T_old[xb1][j])/2.0/dx;
					Txz=(1.0-wxz)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxz*(-3.0*T_old[i][j]+4.0*T_old[xa][j]-T_old[xa1][j])/2.0/dx;
					Tzf=(1.0-wzf)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzf*(3.0*T_old[i][j]-4.0*T_old[i][zb]+T_old[i][zb1])/2.0/dz;
					Tzz=(1.0-wzz)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzz*(-3.0*T_old[i][j]+4.0*T_old[i][za]-T_old[i][za1])/2.0/dz;

					T_xmin=min(T_old[i][j]-dx*Txf,T_old[i][j]+dx*Txz);
					T_zmin=min(T_old[i][j]-dz*Tzf,T_old[i][j]+dz*Tzz);

					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+s[i][j]*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);
				}

			for(i=nx-1;i>=0;i--)
				for(j=0;j<nz;j++)
				{
					xa=i+1;xb=i-1;za=j+1;zb=j-1;
					if(xa==nx) xa=i;
					if(xb==-1) xb=i;
					if(za==nz) za=j;
					if(zb==-1) zb=j;
					
					xa1=i+2;xb1=i-2;za1=j+2;zb1=j-2;
					if(xa1==nx) xa1=xa;
					if(xb1==-1) xb1=xb;
					if(za1==nz) za1=za;
					if(zb1==-1) zb1=zb;
					if(xa1==nx+1) xa1=i;
					if(xb1==-2) xb1=i;
					if(za1==nz+1) za1=j;
					if(zb1==-2) zb1=j;
				
					rxf=(epsi+pow(T_old[i][j]-2*T_old[xb][j]+T_old[xb1][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rxz=(epsi+pow(T_old[xa1][j]-2*T_old[xa][j]+T_old[i][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rzf=(epsi+pow(T_old[i][j]-2*T_old[i][zb]+T_old[i][zb1],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					rzz=(epsi+pow(T_old[i][za1]-2*T_old[i][za]+T_old[i][j],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					wxf=1.0/(1.0+2.0*pow(rxf,2));
					wxz=1.0/(1.0+2.0*pow(rxz,2));
					wzf=1.0/(1.0+2.0*pow(rzf,2));
					wzz=1.0/(1.0+2.0*pow(rzz,2));
					Txf=(1.0-wxf)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxf*(3.0*T_old[i][j]-4.0*T_old[xb][j]+T_old[xb1][j])/2.0/dx;
					Txz=(1.0-wxz)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxz*(-3.0*T_old[i][j]+4.0*T_old[xa][j]-T_old[xa1][j])/2.0/dx;
					Tzf=(1.0-wzf)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzf*(3.0*T_old[i][j]-4.0*T_old[i][zb]+T_old[i][zb1])/2.0/dz;
					Tzz=(1.0-wzz)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzz*(-3.0*T_old[i][j]+4.0*T_old[i][za]-T_old[i][za1])/2.0/dz;

					T_xmin=min(T_old[i][j]-dx*Txf,T_old[i][j]+dx*Txz);
					T_zmin=min(T_old[i][j]-dz*Tzf,T_old[i][j]+dz*Tzz);

					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+s[i][j]*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);
				}
	
			for(i=nx-1;i>=0;i--)
				for(j=nz-1;j>=0;j--)
				{
					xa=i+1;xb=i-1;za=j+1;zb=j-1;
					if(xa==nx) xa=i;
					if(xb==-1) xb=i;
					if(za==nz) za=j;
					if(zb==-1) zb=j;
					
					xa1=i+2;xb1=i-2;za1=j+2;zb1=j-2;
					if(xa1==nx) xa1=xa;
					if(xb1==-1) xb1=xb;
					if(za1==nz) za1=za;
					if(zb1==-1) zb1=zb;
					if(xa1==nx+1) xa1=i;
					if(xb1==-2) xb1=i;
					if(za1==nz+1) za1=j;
					if(zb1==-2) zb1=j;
				
					rxf=(epsi+pow(T_old[i][j]-2*T_old[xb][j]+T_old[xb1][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rxz=(epsi+pow(T_old[xa1][j]-2*T_old[xa][j]+T_old[i][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rzf=(epsi+pow(T_old[i][j]-2*T_old[i][zb]+T_old[i][zb1],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					rzz=(epsi+pow(T_old[i][za1]-2*T_old[i][za]+T_old[i][j],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					wxf=1.0/(1.0+2.0*pow(rxf,2));
					wxz=1.0/(1.0+2.0*pow(rxz,2));
					wzf=1.0/(1.0+2.0*pow(rzf,2));
					wzz=1.0/(1.0+2.0*pow(rzz,2));
					Txf=(1.0-wxf)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxf*(3.0*T_old[i][j]-4.0*T_old[xb][j]+T_old[xb1][j])/2.0/dx;
					Txz=(1.0-wxz)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxz*(-3.0*T_old[i][j]+4.0*T_old[xa][j]-T_old[xa1][j])/2.0/dx;
					Tzf=(1.0-wzf)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzf*(3.0*T_old[i][j]-4.0*T_old[i][zb]+T_old[i][zb1])/2.0/dz;
					Tzz=(1.0-wzz)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzz*(-3.0*T_old[i][j]+4.0*T_old[i][za]-T_old[i][za1])/2.0/dz;

					T_xmin=min(T_old[i][j]-dx*Txf,T_old[i][j]+dx*Txz);
					T_zmin=min(T_old[i][j]-dz*Tzf,T_old[i][j]+dz*Tzz);

					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+s[i][j]*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);
				}


			for(i=0;i<nx;i++)
				for(j=nz-1;j>=0;j--)
				{
					xa=i+1;xb=i-1;za=j+1;zb=j-1;
					if(xa==nx) xa=i;
					if(xb==-1) xb=i;
					if(za==nz) za=j;
					if(zb==-1) zb=j;
					
					xa1=i+2;xb1=i-2;za1=j+2;zb1=j-2;
					if(xa1==nx) xa1=xa;
					if(xb1==-1) xb1=xb;
					if(za1==nz) za1=za;
					if(zb1==-1) zb1=zb;
					if(xa1==nx+1) xa1=i;
					if(xb1==-2) xb1=i;
					if(za1==nz+1) za1=j;
					if(zb1==-2) zb1=j;
				
					rxf=(epsi+pow(T_old[i][j]-2*T_old[xb][j]+T_old[xb1][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rxz=(epsi+pow(T_old[xa1][j]-2*T_old[xa][j]+T_old[i][j],2))/(epsi+pow(T_old[xa][j]-2*T_old[i][j]+T_old[xb][j],2));
					rzf=(epsi+pow(T_old[i][j]-2*T_old[i][zb]+T_old[i][zb1],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					rzz=(epsi+pow(T_old[i][za1]-2*T_old[i][za]+T_old[i][j],2))/(epsi+pow(T_old[i][za]-2*T_old[i][j]+T_old[i][zb],2));
					wxf=1.0/(1.0+2.0*pow(rxf,2));
					wxz=1.0/(1.0+2.0*pow(rxz,2));
					wzf=1.0/(1.0+2.0*pow(rzf,2));
					wzz=1.0/(1.0+2.0*pow(rzz,2));
					Txf=(1.0-wxf)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxf*(3.0*T_old[i][j]-4.0*T_old[xb][j]+T_old[xb1][j])/2.0/dx;
					Txz=(1.0-wxz)*(T_old[xa][j]-T_old[xb][j])/2.0/dx+wxz*(-3.0*T_old[i][j]+4.0*T_old[xa][j]-T_old[xa1][j])/2.0/dx;
					Tzf=(1.0-wzf)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzf*(3.0*T_old[i][j]-4.0*T_old[i][zb]+T_old[i][zb1])/2.0/dz;
					Tzz=(1.0-wzz)*(T_old[i][za]-T_old[i][zb])/2.0/dz+wzz*(-3.0*T_old[i][j]+4.0*T_old[i][za]-T_old[i][za1])/2.0/dz;

					T_xmin=min(T_old[i][j]-dx*Txf,T_old[i][j]+dx*Txz);
					T_zmin=min(T_old[i][j]-dz*Tzf,T_old[i][j]+dz*Tzz);

					if (fabs(T_xmin-T_zmin)>=s[i][j]*dx) T[i][j]=min(T_xmin,T_zmin)+s[i][j]*dx;
					else T[i][j]=(T_xmin+T_zmin+sqrt(2*s[i][j]*dx*s[i][j]*dx-pow(T_xmin-T_zmin,2)))/2;
					T_old[i][j]=min(T[i][j],T_old[i][j]);
				}
			
				lnorm = 0.0;
				for (i=0;i<nx;i++)
				{
					for (j=0;j<nz;j++)
					{							
						lnorm += fabs(T_new[i][j]-T_old[i][j]);							
					}
				}

				for (i=0;i<nx;i++)
				{
					for (j=0;j<nz;j++)
					{							
						T_new[i][j]= T_old[i][j];							
					}
				}									
				loop++;		
					
		}
	}
	end = clock();

	printf("Total %d iterations\n", loop);
	printf("The cost of the run time is %f seconds\n", (double)(end-start)/CLOCKS_PER_SEC);
	
	if(method == 1) sprintf(name,"./traveltime/traveltime_FSM.dat");
	else if(method == 2) sprintf(name,"./traveltime/traveltime_2ndFSM.dat");
	else sprintf(name,"./traveltime/traveltime_stagered_FSM.dat");

	fp=fopen(name,"wb+");
	for(i=0;i<nx;i++)
		for(j=0;j<nz;j++)
		{
			fwrite(&T_old[i][j],sizeof(float),1,fp);
		}
	fclose(fp);

	if(analy == 1)
	{
		sprintf(name, "./traveltime/traveltime_homo_analy.bin");
		fp = fopen(name, "wb+");
		for(i=0;i<nx;i++)
			for(j=0;j<nz;j++)
			{
				T_ana[i][j]=dx*sqrt(pow((i-x0),2)+pow((j-z0),2))*(1.0/v[i][j]);
				fwrite(&T_ana[i][j],sizeof(float),1,fp);
			}
		fclose(fp);			
	}
	else if(analy == 2)
	{

		sprintf(name, "./traveltime/traveltime_gradient_analy.bin");
		fp = fopen(name, "wb+");
		for(i=0;i<nx;i++)
			for(j=0;j<nz;j++)
			{
				float temp = 1+0.5*s[i][j]*s0*(G0*G0)*(pow((i-x0)*dx,2)+pow((j-z0)*dz,2));
				T_ana[i][j] = (1.0/sqrt(G0*G0))*log(temp+sqrt(temp*temp-1));
				fwrite(&T_ana[i][j],sizeof(float),1,fp);
			}	
			fclose(fp);

	}

	for(i=0;i<nx;i++)
	{
		delete [] s[i];delete [] T_old[i];delete [] T_ana[i];delete [] T_new[i];delete [] T[i];delete [] delta[i];delete [] eta[i];delete [] v[i];delete [] vt[i];
	}
	delete [] s;delete [] T_old;delete [] T_ana; delete [] T_new;delete [] T;delete [] delta;delete [] eta;delete [] v;delete [] vt;

}

