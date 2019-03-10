//Adam Trotter, 9 August 2010: Symmetric TRF Statistic
//Requires that both MODEL and DERIVATIVE functions be specified in input file header
//Also requires two slop parameters, named SlopX and SlopY
//
#include <math.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>


#include "../Point.h"
#include "../ExternalResource.h"
#include "nrutil.h"

using namespace std;

typedef double (*MODEL)(vector <double> &, Point & );
typedef double (*CHIMOD)(vector <double> &, Point &, Point &);
MODEL f, dfdx, dfdx2;
int SlopX_Position, SlopY_Position;

StatisticLevelInterface SLI;

extern "C"
{
#ifdef _WIN32
	_declspec(dllexport) int Init(StatisticLevelInterface &S)
#else
	int Init(StatisticLevelInterface &S)
#endif
	{
		SLI = S;
		ExternalResource E;
		E.SetSLI(S);
		E.Load(SLI.Statistic_Arguments["MODEL"].at(0).c_str(), SLI.Statistic_Arguments["MODEL"].at(1).c_str());
		f = (MODEL) E.Address;
		printf("Model Loaded\n");
		E.Load(SLI.Statistic_Arguments["DERIVATIVE"].at(0).c_str(), SLI.Statistic_Arguments["DERIVATIVE"].at(1).c_str());
		dfdx = (MODEL) E.Address;
		printf("First Derivative Loaded\n");
		E.Load(SLI.Statistic_Arguments["DERIVATIVE2"].at(0).c_str(), SLI.Statistic_Arguments["DERIVATIVE2"].at(1).c_str());
		dfdx2 = (MODEL) E.Address;
		printf("Second Derivative Loaded\n");

		SlopX_Position = SLI.Parameter_Positions["SlopX"];
		SlopY_Position = SLI.Parameter_Positions["SlopY"];



		return 0;
	}



#ifdef _WIN32
	_declspec(dllexport) double TRF(vector<vector<double> > &Organism)
#else
	double TRF(vector<vector<double> > &Organism)
#endif
	{
		double slopx, slopy;
		double xn,sigxn,yn,sigyn,weight;
		double ex,ey,e0,arg,fac1;
		double chisqrold,chisqr,chisqrmin,tol,xt,xtold,yt,yxt,ytold,xtmin,ytmin,slope,slopemin;
		int zount,j,n;
		double chi = 0.0;
		Point* p;
		Point pt;

		tol = 1.0E-8;

		for(j=0;j<(signed)SLI.Data.size();j++)
		{
			slopx = Organism[j][SlopX_Position]; 
			slopy = Organism[j][SlopY_Position];
			for(n=0;n<(signed)SLI.Data[j].size();n++)
			{
				p = &SLI.Data[j][n];
				pt = *p;
				//pt.X[0] = p->X[0];
				//pt.Y = p->Y;

				//pt = *p;
				xn = p->X[0]; // data point x-value
				sigxn = p->X_PLUS[0]; // data point 1-sigma x-errorbar (assume symmetric Gaussian)
				yn = p->Y; // data point y-value
				sigyn = p->Y_PLUS; // data point 1-sigma y-errorbar (assume symmetric Gaussian)
				weight = p->X[1]; // data point weight (should be in second "X column", after data X errorbar
				ex = sqrt(sigxn*sigxn+slopx*slopx);
				ey = sqrt(sigyn*sigyn+slopy*slopy);
				
				//Determine tangent point (xt,yt) of model wrt convolved error ellipse
				chisqrold = 1;
				chisqr = 2.0E6;
				chisqrmin = 1.0E6;
				zount = 0;
				xt = xn;
				yt = yn;
				
				while(sqrt((chisqr-chisqrold)*(chisqr-chisqrold)) > tol && zount < 100)
				{
					zount++;
					yxt = (f)(Organism[j],pt);
					slope = (dfdx)(Organism[j],pt);
					
					xtold = xt;
					if (slope == 0.0)
						xt = xn;
					else
						xt = (slope*xt + ey*ey/ex/ex/slope*xn + yn - yxt) / (slope+ey*ey/ex/ex/slope);
					ytold = yt;				
					yt = yxt + slope * (xt - xtold);
					pt.X[0] = xt;
					pt.Y = yt;
					chisqrold = chisqr;
					e0 = sqrt(slope*slope*ex*ex + ey*ey);
					arg = yn-(yt+slope*(xn-xt));
					//arg = -0.5*arg*arg/e0/e0;
					//chisqr = -2.0*arg;
					chisqr = arg*arg/e0/e0;
					if (chisqr < chisqrmin) //in case tangent point iteration fails to converge, keep the minimum as a backup
					{
						chisqrmin = chisqr;
						xtmin = xt;
						ytmin = yt;
						slopemin = slope;
					}
				}

				chisqr = chisqrmin;
				slope = slopemin;

				//TRF Prefactor
				if (ex == 0)
					fac1 = 1.0;
				else
					fac1 = (slope*slope+ey*ey/ex/ex)/sqrt(slope*slope+pow(ey/ex,4));
				
				// -2ln(prob)
				chi += weight*(chisqr - 2.0*log(fac1) + 2.0*log(e0)); 
			}

		}
		//Return -2*ln(prob).  Prior functions do the same.
		return chi;
	}

#ifdef _WIN32
	_declspec(dllexport) double TRF_Line(vector<vector<double> > &Organism)
#else
	double TRF_Line(vector<vector<double> > &Organism)
#endif
	{
		double slopx, slopy, scx, scy, sc, fx, fy;
		double xn,sigxn,yn,sigyn,weight;
		double ex,ey,e0,arg,fac1;
		double chisqr,xt,yt,slope;
		double b,m,xp;
		int j,n;
		double chi = 0.0;
		Point *p;
		//Point pt;
		/*
		double sumwt = 0.0;
		double sumlnsig = 0.0;
		double sumterm = 0.0;
		double sumlnfac = 0.0;
		double sumfac = 0.0;
		double sumchi = 0.0;
		double sumsig = 0.0;
		double sumr = 0.0;
		double sumr2 = 0.0;
		double sumr4 = 0.0;
		double summ2 = 0.0;
		double sumar = 0.0;
		double sumbr = 0.0;
		double sumax = 0.0;
		double sumbx = 0.0;
		double sumay = 0.0;
		double sumby = 0.0;
		double r,ai,bi;
		*/

		for(j=0;j<(signed)SLI.Data.size();j++)
		{
			//Linear Model: Organism c[0]=b, c[1]=theta, c[2]=pivot poitn, c[3]=SlopX, c[4]=SlopY, c[5]=ScaleX, c[6]=ScaleY
			slopx = Organism[j][SlopX_Position]; 
			slopy = Organism[j][SlopY_Position];
			//slopx = pow(10.0,Organism[j][SlopX_Position]);
			//slopy = pow(10.0,Organism[j][SlopY_Position]);
			b = Organism[j][0];
			m = tan(Organism[j][1]*0.0174532925);
			//m = Organism[j][1];
			xp = Organism[j][2];
			//scx = Organism[j][5];
			//scy = Organism[j][6];
			scx = pow(10.0,Organism[j][5]);
			scy = pow(10.0,Organism[j][6]);
			fx = pow(10.0,Organism[j][7]);
			fy = pow(10.0,Organism[j][8]);
			sc = scy/scx;
			//cout << " b = " << b << " m = " << m << endl;
			for(n=0;n<(signed)SLI.Data[j].size();n++)
			{
				p = &SLI.Data[j][n];
				//pt = *p;
				xn = p->X[0]; // data point x-value
				sigxn = p->X_PLUS[0]; // data point 1-sigma x-errorbar (assume symmetric Gaussian)
				yn = p->Y; // data point y-value
				sigyn = p->Y_PLUS; // data point 1-sigma y-errorbar (assume symmetric Gaussian)
				sigxn = sigxn * fx;
				sigyn = sigyn * fy;
				//cout << "xn = " << xn << " yn = " << yn << endl;
				//yn = p->X[2];
				//sigyn = p->X_PLUS[2]; //Kluge to get y-errorbars in...they're not being read from the final column in the input file!
				

				// Option 1: Scale data points and intrinsic errorbars by scx, scy
				//xn = xn*scx;
				//sigxn = sigxn*scx;
				//yn = yn*scy;
				//sigyn = sigyn*scy;
				// Option 2: Don't scale data, but include sc=scy/scx in prefactor (below)
				weight = p->X[1]; // data point weight (should be in second "X column", after data X errorbar
				//cout << "xn = " << xn << " sigxn = " << sigxn << " yn = " << yn << " sigyn = " << sigyn <<  " weight = " << weight << endl;
				//system("pause");

				ex = sqrt(sigxn*sigxn+slopx*slopx);
				ey = sqrt(sigyn*sigyn+slopy*slopy);
				
				//Determine tangent point (xt,yt) of model wrt convolved error ellipse
				xt = (m*ex*ex*(yn-b+m*xp)+ey*ey*xn)/(m*m*ex*ex+ey*ey);
				yt = b+m*(xt-xp);
				slope = m;
				e0 = sqrt(slope*slope*ex*ex + ey*ey);
				arg = yn-(yt+slope*(xn-xt));
				//arg = yn-slope*(xn-xp)-b;
				chisqr = arg*arg/e0/e0;

				//TRF Prefactor
				if (ex == 0)
					fac1 = 1.0;
				else
					//fac1 = (slope*slope+ey*ey/ex/ex)/sqrt(slope*slope+sc*sc*pow(ey/ex,4))/scx;
					fac1 = (slope*slope+ey*ey/ex/ex)/sqrt(slope*slope+sc*sc*pow(ey/ex,4));
				
				// -2ln(prob)
				// -2ln(prob)
				//DO5 Prefactor
				//fac1 = 1.0;
				
				chi += weight*(chisqr - 2.0*log(fac1) + 2.0*log(e0));
				
				/*
				r = ey/ex;
				ai = 4.0*sc*sc*r*r*r/(slope*slope+sc*sc*r*r*r*r)-4.0*r/(slope*slope+r*r);
				bi = 1.0/e0/e0-chisqr/e0/e0;
				sumwt += weight;
				sumlnsig += weight*log(e0);
				sumsig += weight*e0;
				sumr += weight*r;
				sumr2 += weight*r*r;
				sumr4 += weight*r*r*r*r;
				summ2 += weight*slope*slope;
				sumchi += weight*chisqr;
				sumlnfac += weight*log((slope*slope+r*r)/sqrt(slope*slope+sc*sc*r*r*r*r));
				sumfac += weight*(slope*slope+r*r)/sqrt(slope*slope+sc*sc*r*r*r*r);
				sumar += weight*ai;
				sumbr += weight*2.0*slopx*slopy*(1.0-slope*slope/r/r)*bi;
				sumax += -weight*slopy/slopx/slopx*ai;
				sumbx += weight*2.0*slope*slope*slopx*bi;
				sumay += weight/slopx*ai;
				sumby += weight*2.0*slopy*bi;
				*/
			}

		}
		//Return -2*ln(prob). Sign gets changed in BayesianInference.dll
		/*
		cout << "chi = " << chi << endl;
		cout << " N = " << sumwt << endl;
		cout << " SUMlnfac = " << sumlnfac << endl;
		cout << " SUMlnSIG = " << sumlnsig << endl;
		cout << " SUMchi = " << sumchi << endl;
		cout << " <r> = " << sumr/sumwt << endl;
		cout << " <r^2> = " << sumr2/sumwt << endl;
		cout << " <r^4> = " << sumr4/sumwt << endl;
		cout << " <m^2> = " << summ2/sumwt << endl;
		cout << " <SIG> = " << sumsig/sumwt << endl;
		cout << " <fac> = " << sumfac/sumwt << endl;
		cout << " <Ar> = " << sumar/sumwt << endl;
		cout << " <Br> = " << sumbr/sumwt << endl;
		cout << " <Ax> = " << sumax/sumwt << endl;
		cout << " <Bx> = " << sumbx/sumwt << endl;
		cout << " <Ay> = " << sumay/sumwt << endl;
		cout << " <By> = " << sumby/sumwt << endl;
		system("pause");
		*/
		return chi;

	}


#ifdef _WIN32
	_declspec(dllexport) void mnbrak(CHIMOD func, vector<double> c, Point pn, double &ax, double &bx, double &cx)
#else
	void mnbrak(CHIMOD func, vector<double> c, Point pn, double &ax, double &bx, double &cx)
#endif
		{
			double GOLD = 1.618034;
			double GLIMIT = 100.0;
			double TINY = 1.0E-20;
			double ulim,u,r,q,fu,dum,fu2;
			double fa,fb,fc,ya,yb,yc,yu;
			Point pt;

			pt.X.push_back(0.0);
			pt.X[0]=ax;
			//ya = (f)(c,pt);
			//fa = (ax-xn)*(ax-xn)/ex/ex+(ya-yn)*(ya-yn)/ey/ey;
			fa = (func)(c, pt, pn);
			pt.X[0]=bx;
			//yb = (f)(c,pt);
			//fb = (bx-xn)*(bx-xn)/ex/ex+(yb-yn)*(yb-yn)/ey/ey;
			fb = (func)(c, pt, pn);
			//cout << "mnbrak: ax = " << ax << " fa = " << fa << " bx = " << bx << " fb = " << fb << endl;
			if (fb > fa)
			{
				SHFT(dum,ax,bx,dum)
				SHFT(dum,fb,fa,dum)
			}
			cx = bx+GOLD*(bx-ax);
			pt.X[0]=cx;
			//yc = (f)(c,pt);
			//fc = (cx-xn)*(cx-xn)/ex/ex+(yc-yn)*(yc-yn)/ey/ey;
			fc = (func)(c, pt, pn);
			//cout << "mnbrak: cx = " << cx << " fc = " << fc << endl;
			while (fb > fc)
			{
				r=(bx-ax)*(fb-fc);
				q=(bx-cx)*(fb-fa);
				u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0+SIGN(DMAX(fabs(q-r),TINY),q-r));
				ulim=bx+GLIMIT*(cx-bx);
				if ((bx-u)*(u-cx) > 0.0)
				{
					pt.X[0]=u;
					//yu = (f)(c,pt);
					//fu = (u-xn)*(u-xn)/ex/ex+(yu-yn)*(yu-yn)/ey/ey;
					fu = (func)(c, pt, pn);
					if (fu < fc)
					{
						ax=bx;
						bx=u;
						fa=fb;
						fb=fu;
						//cout << "mnbrak: ax = " << ax << " bx = " << bx << " cx = " << cx << endl;
						//cout << "mnbrak: fa = " << fa << " fb = " << fb << " fc = " << fc << endl;
						return;
					}
					else if (fu > fb)
					{
						cx=u;
						fc=fu;
						//cout << "mnbrak: ax = " << ax << " bx = " << bx << " cx = " << cx << endl;
						//cout << "mnbrak: fa = " << fa << " fb = " << fb << " fc = " << fc << endl;
						return;
					}
					u=cx+GOLD*(cx-bx);
					pt.X[0]=u;
					//yu = (f)(c,pt);
					//fu = (u-xn)*(u-xn)/ex/ex+(yu-yn)*(yu-yn)/ey/ey;
					fu = (func)(c, pt, pn);
				}
				else if ((cx-u)*(u-ulim) > 0.0)
				{
					pt.X[0]=u;
					//yu = (f)(c,pt);
					//fu = (u-xn)*(u-xn)/ex/ex+(yu-yn)*(yu-yn)/ey/ey;
					fu = (func)(c, pt, pn);
					if (fu < fc)
					{
						SHFT(bx,cx,u,cx+GOLD*(cx-bx))
						pt.X[0]=u;
						//yu = (f)(c,pt);
						//fu2 = (u-xn)*(u-xn)/ex/ex+(yu-yn)*(yu-yn)/ey/ey;
						fu2 = (func)(c, pt, pn);
						SHFT(fb,fc,fu,fu2);
					}
				}
				else if ((u-ulim)*(ulim-cx) >= 0.0)
				{
					u=ulim;
					pt.X[0]=u;
					//yu = (f)(c,pt);
					//fu = (u-xn)*(u-xn)/ex/ex+(yu-yn)*(yu-yn)/ey/ey;
					fu = (func)(c, pt, pn);
				}
				else
				{
					u=cx+GOLD*(cx-bx);
					pt.X[0]=u;
					//yu = (f)(c,pt);
					//fu = (u-xn)*(u-xn)/ex/ex+(yu-yn)*(yu-yn)/ey/ey;
					fu = (func)(c, pt, pn);
				}
				SHFT(ax,bx,cx,u)
				SHFT(fa,fb,fc,fu)
			}			

		}

#ifdef _WIN32
	_declspec(dllexport) double brent(CHIMOD func, vector<double> c, Point pn, double ax, double bx, double cx)
#else
	double brent(CHIMOD func, vector<double> c, Point pn, double ax, double bx, double cx)
#endif
		{
			int ITMAX = 100;
			double CGOLD = 0.3819660;
			double ZEPS = 1.0E-10;
			double tol = 1.0E-8;
			int iter;
			double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
			double e=0.0;
			double yx,yu,xmin;
			Point pt;

			pt.X.push_back(0.0);
			
			a=(ax < cx ? ax : cx);
			b=(ax > cx ? ax : cx);
			x=w=v=bx;
			pt.X[0]=x;
			//yx = (f)(c,pt);
			//fx = (x-xn)*(x-xn)/ex/ex+(yx-yn)*(yx-yn)/ey/ey;
			fx = (func)(c, pt, pn);
			fw=fv=fx;
			for (iter=1;iter<=ITMAX;iter++)
			{
				//cout << "iter = " << iter << endl;
				xm=0.5*(a+b);
				tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
				if (fabs(x-xm) <= (tol2-0.5*(b-a)))
				{
					xmin=x;
					return xmin;
				}
				if (fabs(e) > tol1)
				{
					r=(x-w)*(fx-fv);
					q=(x-v)*(fx-fw);
					p=(x-v)*q-(x-w)*r;
					q=2.0*(q-r);
					if (q > 0.0) p = -p;
					q=fabs(q);
					etemp=e;
					e=d;
					if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
						d=CGOLD*(e=(x >= xm ? a-x : b-x));
					else
					{
						d=p/q;
						u=x+d;
						if (u-a < tol2 || b-u < tol2)
							d=SIGN(tol1,xm-x);
					}
				}
				else
				{
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				}
				u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
				pt.X[0]=u;
				//yu = (f)(c,pt);
				//fu = (u-xn)*(u-xn)/ex/ex+(yu-yn)*(yu-yn)/ey/ey;
				fu = (func)(c, pt, pn);
				if (fu <= fx)
				{
					if(u >= x) a=x; else b=x;
					SHFT(v,w,x,u);
					SHFT(fv,fw,fx,fu);
				}
				else
				{
					if (u < x) a=u; else b=u;
					if (fu <= fw || w == x)
					{
						v=w;
						w=u;
						fv=fw;
						fw=fu;
					}
					else if (fu <= fv || v == x || v == w)
					{
						v=u;
						fv=fu;
					}
				}
			}
			//cout << "Too many iterations in brent" << endl;
			xmin = -3.0;
			return xmin;
		}

#ifdef _WIN32
	_declspec(dllexport) double rtbis(CHIMOD func, vector<double> c, Point pn, double x1, double x2, double xacc)
#else
	double rtbis(CHIMOD func, vector<double> c, Point pn, double x1, double x2, double xacc)
#endif
		{
			int j;
			int JMAX=100;
			double dx,f,fmid,xmid,rtb;
			Point pt;
			pt.X.push_back(0.0);
			pt.X[0]=x1;
			f=(func)(c, pt, pn);
			pt.X[0]=x2;
			fmid=(func)(c, pt, pn);
			/*
			if (f*fmid >= 0.0) 
			{
				
				cout << "Root must be bracketed in rtbis" << endl;
				cout << c[0] << endl;
				cout << c[1] << endl;
				cout << c[2] << endl;
				cout << c[3] << endl;
				cout << c[4] << endl;
				cout << c[5] << endl;
				cout << c[6] << endl;
				cout << c[7] << endl;
				cout << c[8] << endl;
				cout << c[9] << endl;
				cout << c[10] << endl;
				cout << c[11] << endl;
				cout << c[12] << endl;
				cout << " x1 = " << x1 << " f = " << f << " x2 = " << x2 << " fmid = " << fmid << endl;
				system("pause");
				
			}
			*/
			rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
			for (j=1;j<JMAX;j++)
			{
				xmid=rtb+(dx *= 0.5);
				pt.X[0]=xmid;
				fmid = (func)(c, pt, pn);
				if (fmid <= 0.0) rtb=xmid;
				if (fabs(dx) < xacc || fmid == 0.0) return rtb;
			}
			//cout << "Too many bisections in rtbis!" << endl;
			return 0.0;
	}


#ifdef _WIN32
	_declspec(dllexport) double DCHI(vector<double> &c, Point & pt, Point & pn)
#else
	double DCHI(vector<double> &c, Point & pt, Point & pn)
#endif
		{
			//DCHI = (xt-xn)^2/ex^2+(yt-yn)^2/ey^2.  Here, multiplied by (ex*ey)^2 to avoid divide by zero
			double yt = (f)(c, pt);
			return (pt.X[0]-pn.X[0])*(pt.X[0]-pn.X[0])*pn.Y_PLUS*pn.Y_PLUS+(yt-pn.Y)*(yt-pn.Y)*pn.X_PLUS[0]*pn.X_PLUS[0];
		}

#ifdef _WIN32
	_declspec(dllexport) double DCHIDX(vector<double> &c, Point & pt, Point & pn)
#else
	double DCHIDX(vector<double> &c, Point &pt, Point &pn)
#endif
		{
			//DCHIDX = d/dx(DCHI).  Multiplied by (ex*ey)^2 to avoid divide by zero
			double yt = (f)(c, pt);
			double dydx = (dfdx)(c, pt);
			return (pt.X[0]-pn.X[0])*pn.Y_PLUS*pn.Y_PLUS+(yt-pn.Y)*dydx*pn.X_PLUS[0]*pn.X_PLUS[0];
		}

#ifdef _WIN32
	_declspec(dllexport) double DCHIDX2(vector<double> &c, Point & pt, Point & pn)
#else
	double DCHIDX2(vector<double> &c, Point &pt, Point &pn)
#endif
		{
			//DCHIDX2 = second derivative d2/dx2(DCHI).  Multiplied by (ex*ey)^2 to avoid divide by zero
			double yt = (f)(c, pt);
			double dydx = (dfdx)(c, pt);
			double dydx2 = (dfdx2)(c, pt);
			return ((yt-pn.Y)*dydx2+dydx*dydx)*pn.X_PLUS[0]*pn.X_PLUS[0]+pn.Y_PLUS*pn.Y_PLUS;
		}

#ifdef _WIN32
	_declspec(dllexport) double DXY(vector<double> &c, Point & pt, Point & pn)
#else
	double DXY(vector<double> &c, Point &pt, Point &pn)
#endif
		{
			return pn.Y-(f)(c,pt);
		}

#ifdef _WIN32
	_declspec(dllexport) double TRF_BHc2(vector<vector<double> > &Organism)
#else
	double TRF_BHc2(vector<vector<double> > &Organism)
#endif
	{
		//c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s, 
		//c[7]=SlopX, c[8]=Slopy, c[9]=logScaleX, c[10]=logScaleY, c[11]=logfX, c[12]=logfY;

		double slopx, slopy;
		double xn,sigxn,yn,sigyn,weight;
		double ex,ey,e0,arg,fac;
		double chisqr,xt,yt,yxt,slope;
		double xint,yint,b1,m1,p1,b2,m2,p2,s;
		int j,n;
		double chi = 0.0;
		double ax,bx,cx,xmin;
		double f1,f2,x0,y0,x1,x2,yt1,yt2,arg1,arg2,xt1,xt2,slope1,slope2,chisqr1,chisqr2,yxn;
		double e01,e02,fac1,fac2;
		double scx,scy,sc,fx,fy;
		double dx,fmid,f0,xmid,xacc,xdchi2min,dchi2min;
		Point* p;
		Point pt, px1, px2, pn;
		CHIMOD dchi,dchidx,dchidx2;
		dchi = &DCHI;
		dchidx = &DCHIDX;
		dchidx2 = &DCHIDX2;
		pt.X.push_back(0.0);
		xacc = 1.0E-6;
	
		/*
		double sumwt = 0.0;
		double sumlnsig = 0.0;
		double sumterm = 0.0;
		double sumlnfac = 0.0;
		double sumfac = 0.0;
		double sumchi = 0.0;
		double sumsig = 0.0;
		double sumr = 0.0;
		double sumr2 = 0.0;
		double sumr4 = 0.0;
		double summ2 = 0.0;
		double sumar = 0.0;
		double sumbr = 0.0;
		double sumax = 0.0;
		double sumbx = 0.0;
		double sumay = 0.0;
		double sumby = 0.0;
		double r,ai,bi;
		*/
		//string filename;
		//char* file;

		//filename = "diagnose.txt";
		//file = const_cast<char*> (filename.c_str());
		//std::ofstream tmp;
		//tmp.open(name);
		//ofstream outfile("diagnose.txt");
		for(j=0;j<(signed)SLI.Data.size();j++)
		{
			slopx = Organism[j][SlopX_Position]; 
			slopy = Organism[j][SlopY_Position];
			//intersection of 2 lines: c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s
			b1 = Organism[j][0];
			m1 = tan(Organism[j][1]*0.0174532925);
			p1 = Organism[j][2];
			b2 = Organism[j][3];
			m2 = tan(Organism[j][4]*0.0174532925);
			p2 = Organism[j][5];
			s = Organism[j][6];
			//c = Organism[j];
			scx = pow(10.0,Organism[j][9]); //x-scale factor for data points
			scy = pow(10.0,Organism[j][10]); //y-scale factor for data points
			sc = scy/scx;
			fx = pow(10.0,Organism[j][11]); //Scaling factor for x-errorbars.  fx<=1
			fy = pow(10.0,Organism[j][12]); //  and for y-errorbars
			//Intercept of two lines b1+m1(xint-p1) = b2+m2(xint-p2)
			xint = (b2-m2*p2-b1+m1*p1)/(m1-m2);
			//For BHc2, or RVc2 with theta2>0, need to find where derivative of model=0, to help bracket tangent points  
			if (m1*m2 >= 0) // if slopes both have the same sign, there is no place where dydx=0; use intersection point
							//  Note: this should never happen for BH vs. c2!
			{
				x0 = xint; 
			}
			else // x0 = point where dydx=0
			{
				if(m1<0)
				{
					x0 = (log(-m1)-log(m2)-s*b2+s*m2*p2+s*b1-s*m1*p1)/(s*(m2-m1));
				}
				else
				{
					x0 = (log(-m2)-log(m1)-s*b1+s*m1*p1+s*b2-s*m2*p2)/(s*(m1-m2));
				}
			}
			pt.X[0]=x0;
			y0 = (f)(Organism[j],pt); // y-value of peak of BH vs. c2 curve
			
			for(n=0;n<(signed)SLI.Data[j].size();n++)
			{
				p = &SLI.Data[j][n];
				pt = *p;
				pn = pt;
				xn = p->X[0]; // data point x-value
				sigxn = p->X_PLUS[0]; // data point 1-sigma x-errorbar (assume symmetric Gaussian)
				yn = p->Y; // data point y-value
				sigyn = p->Y_PLUS; // data point 1-sigma y-errorbar (assume symmetric Gaussian)
				sigxn = sigxn*fx;
				sigyn = sigyn*fy;
				weight = p->X[1]; // data point weight (should be in second "X column", after data X errorbar
				//cout << "xn = " << xn << " sigxn = " << sigxn << " yn = " << yn << " sigyn = " << sigyn << " weight = " << weight << endl;
				//system("pause");
				ex = sqrt(sigxn*sigxn+slopx*slopx);
				ey = sqrt(sigyn*sigyn+slopy*slopy);
				//pn = Data Point with errorbars convolved with slop.
				pn.X_PLUS[0]=ex;
				pn.Y_PLUS=ey; 
				yxn = (f)(Organism[j],pn); //Value of model curve at x=xn.
				//cout << "xn = " << xn << " ex = " << ex << " yn = " << yn << " ey = " << ey << endl;
				//cout << "x0 = " << x0 << " y0 = " << y0 << " yxn = " << yxn << endl;

				if (yn > y0) // If data point is above peak of curve, there is only one tangent point, on the same side of the
							 //  peak as the data point...find it using rtbis on DCHIDX
				{
					//cout << " yn > y0: " << endl;
					if (xn<x0)
					{
						xt = rtbis(dchidx,Organism[j],pn,-1.0,3.0,xacc);
					}
					else
					{
						xt = rtbis(dchidx,Organism[j],pn,-1.0,3.0,xacc);
					}
					pt.X[0]=xt;
					yt = (f)(Organism[j],pt);
					slope = (dfdx)(Organism[j],pt);
					//cout << " xt = " << xt << " yt = " << yt << " slope = " << slope << endl;
				}

				else // If data point is below peak of curve, there can be up to 3 tangent points
								   // In any case, 2nd derivative DCHIDX2 will have a single minimum, somewhere near x0.
								   //  Find it using mnbrak and brent
								   // If the minimum is > 0, there is only one tangent point
								   // Otherwise, there may be more than one
				{
					//cout << " yn < y0:" << endl;
					ax = x0;
					bx = x0+0.1;
					mnbrak(dchidx2,Organism[j],pn,ax,bx,cx);
					xdchi2min = brent(dchidx2,Organism[j],pn,ax,bx,cx);
					pt.X[0]=xdchi2min;
					dchi2min = (dchidx2)(Organism[j],pt,pn);
					//cout << " xdchi2min = " << xdchi2min << " dchi2min = " << dchi2min << endl;
					if (dchi2min >= 0.0) // If dchi2min > 0, there is only one tangent point, on the same side of the peak as the data point
						                 //   Find it using rtbis on dchidx
					{
						if(xn > x0)
						{
							xt = rtbis(dchidx, Organism[j], pn, -1.0, 3.0, xacc);
						}
						else
						{
							xt = rtbis(dchidx, Organism[j], pn, -1.0, 3.0, xacc);
						}
						pt.X[0]=xt;
						yt = (f)(Organism[j],pt);
						slope = (dfdx)(Organism[j],pt);
						//cout << " dchi2min > 0:" << endl;
						//cout << " xt = " << xt << " yt = " << yt << " slope = " << slope << endl;
					}
					else  // If dchi2min < 0, there may be up to three tangent points, one of which is at a local maximimum of dchi
						  // The tangent points, if they exist, are bracketed by zeroes of dchidx2
						  //  Find the bracketing points x1, x2 using rtbis on dchidx2
					{
						//cout << " dchi2min < 0" << endl;
						x1 = rtbis(dchidx2, Organism[j], pn, -1.0, xdchi2min, xacc);
						x2 = rtbis(dchidx2, Organism[j], pn, xdchi2min, 3.0, xacc);
						//cout << " x1 = " << x1 << " x2 = " << x2 << endl;
						// The local maximum of dchi is between x1 and x2...I think.  Ignore it.
						// Find the local minima tangent points < x1 and > x2, if they exist, using rtbis on dchidx
						pt.X[0]=-1.0;
						f1 = (dchidx)(Organism[j],pt,pn);
						pt.X[0]=x1;
						f2 = (dchidx)(Organism[j],pt,pn);
						//cout <<  "f1 = " << f1 << " f2 = " << f2 << endl;
						if (f1*f2 < 0.0)
							xt1 = rtbis(dchidx, Organism[j], pn, -3.0, x1, xacc);
						else
							xt1 = -1.0;
						pt.X[0]=xt1;
						yt1 = (f)(Organism[j],pt);
						slope1 = (dfdx)(Organism[j],pt);
						e01 = sqrt(slope1*slope1*ex*ex + ey*ey);
						arg1 = yn-(yt1+slope1*(xn-xt1));
						if(ex == 0.0)
							fac1 = 1.0;
						else
							fac1 = (slope1*slope1+ey*ey/ex/ex)/sqrt(slope1*slope1+sc*sc*pow(ey/ex,4));
						chisqr1 = arg1*arg1/e01/e01 - 2.0*log(fac1) + 2.0*log(e01);
						
						pt.X[0]=x2;
						f1=(dchidx)(Organism[j],pt,pn);
						pt.X[0]=3.0;
						f2=(dchidx)(Organism[j],pt,pn);
						//cout << " f1 = " << f1 << " f2 = " << f2 << endl;
						if(f1*f2 < 0.0)
							xt2 = rtbis(dchidx, Organism[j], pn, x2, 5.0, xacc);
						else
							xt2 = 3.0;
						pt.X[0]=xt2;
						yt2 = (f)(Organism[j],pt);
						slope2 = (dfdx)(Organism[j],pt);
						e02 = sqrt(slope2*slope2*ex*ex + ey*ey);
						arg2 = yn-(yt2+slope2*(xn-xt2));
						if(ex == 0.0)
							fac2 = 1.0;
						else
							fac2 = (slope2*slope2+ey*ey/ex/ex)/sqrt(slope2*slope2+sc*sc*pow(ey/ex,4));
						chisqr2 = arg2*arg2/e02/e02 - 2.0*log(fac2) + 2.0*log(e02);
						//cout << "xt1 = " << xt1 << " yt1 = " << yt1 << " chisqr1 = " << chisqr1 << endl;
						//cout << "xt2 = " << xt2 << " yt2 = " << yt2 << " chisqr2 = " << chisqr2 << endl;

						//Choose the tangent point that minimizes chisqr
						if(chisqr1<chisqr2)
						{
							xt=xt1;
							yt=yt1;
							slope = slope1;
						}
						else
						{
							xt=xt2;
							yt=yt2;
							slope = slope2;
						}
					}
				}

				e0 = sqrt(slope*slope*ex*ex + ey*ey);
				arg = yn-(yt+slope*(xn-xt));
				if(ex == 0.0)
					fac = 1.0;
				else
					fac = (slope*slope+ey*ey/ex/ex)/sqrt(slope*slope+sc*sc*pow(ey/ex,4));
				chisqr = arg*arg/e0/e0 - 2.0*log(fac) + 2.0*log(e0);
				
				//outfile << xn << " " << ex << " " << yn << " " << ey << " " << xt << " " << yt << " " << weight << " " << slope << " " << chisqr << endl;
				//if(yn>7)
				//	cout << xn << " " << ex << " " << yn << " " << ey << " " << xt << " " << yt << " " << slope << " " << chisqr << endl;
				chi += weight*chisqr;
				//system("pause");
				/*
				r = ey/ex;
				ai = 4.0*sc*sc*r*r*r/(slope*slope+sc*sc*r*r*r*r)-4.0*r/(slope*slope+r*r);
				bi = 1.0/e0/e0-chisqr/e0/e0;
				sumwt += weight;
				sumlnsig += weight*log(e0);
				sumsig += weight*e0;
				sumr += weight*r;
				sumr2 += weight*r*r;
				sumr4 += weight*r*r*r*r;
				summ2 += weight*slope*slope;
				sumchi += weight*chisqr;
				sumlnfac += weight*log((slope*slope+r*r)/sqrt(slope*slope+sc*sc*r*r*r*r));
				sumfac += weight*(slope*slope+r*r)/sqrt(slope*slope+sc*sc*r*r*r*r);
				sumar += weight*ai;
				sumbr += weight*2.0*slopx*slopy*(1.0-slope*slope/r/r)*bi;
				sumax += -weight*slopy/slopx/slopx*ai;
				sumbx += weight*2.0*slope*slope*slopx*bi;
				sumay += weight/slopx*ai;
				sumby += weight*2.0*slopy*bi;
				*/
			}

		}
		//Return -2*ln(prob). Sign gets changed in BayesianInference.dll
		/*
		cout << "chi = " << chi << endl;
		cout << " N = " << sumwt << endl;
		cout << " SUMlnfac = " << sumlnfac << endl;
		cout << " SUMlnSIG = " << sumlnsig << endl;
		cout << " SUMchi = " << sumchi << endl;
		cout << " <r> = " << sumr/sumwt << endl;
		cout << " <r^2> = " << sumr2/sumwt << endl;
		cout << " <r^4> = " << sumr4/sumwt << endl;
		cout << " <m^2> = " << summ2/sumwt << endl;
		cout << " <SIG> = " << sumsig/sumwt << endl;
		cout << " <fac> = " << sumfac/sumwt << endl;
		cout << " <Ar> = " << sumar/sumwt << endl;
		cout << " <Br> = " << sumbr/sumwt << endl;
		cout << " <Ax> = " << sumax/sumwt << endl;
		cout << " <Bx> = " << sumbx/sumwt << endl;
		cout << " <Ay> = " << sumay/sumwt << endl;
		cout << " <By> = " << sumby/sumwt << endl;
		system("pause");
		*/
		//outfile.close();
		//system("pause");
		//cout << chi << endl;
		return chi;
	}

#ifdef _WIN32
	_declspec(dllexport) double TRF_BHc2_X(vector<vector<double> > &Organism)
#else
	double TRF_BHc2_X(vector<vector<double> > &Organism)
#endif
	{
		//c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s, 
		//c[7]=SlopX, c[8]=Slopy, c[9]=logScaleX, c[10]=logScaleY, c[11]=logfX, c[12]=logfY;

		double slopx, slopy;
		double xn,sigxn,yn,sigyn,weight;
		double ex,ey,e0,arg,fac;
		double chisqr,xt,yt,yxt,slope;
		double xint,yint,b1,m1,p1,b2,m2,p2,s;
		int j,n;
		double chi = 0.0;
		double ax,bx,cx,xmin;
		double f1,f2,x0,y0,x1,x2,yt1,yt2,arg1,arg2,xt1,xt2,slope1,slope2,chisqr1,chisqr2,yxn;
		double e01,e02,fac1,fac2;
		double scx,scy,sc,fx,fy;
		double dx,fmid,f0,xmid,xacc,xdchi2min,dchi2min;
		Point* p;
		Point pt, px1, px2, pn;
		CHIMOD dchi,dchidx,dchidx2,dxy;
		dchi = &DCHI;
		dchidx = &DCHIDX;
		dchidx2 = &DCHIDX2;
		dxy = &DXY;
		pt.X.push_back(0.0);
		xacc = 1.0E-12;
	
		for(j=0;j<(signed)SLI.Data.size();j++)
		{
			slopx = Organism[j][SlopX_Position]; 
			slopy = Organism[j][SlopY_Position];
			//intersection of 2 lines: c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s
			b1 = Organism[j][0];
			m1 = tan(Organism[j][1]*0.0174532925);
			p1 = Organism[j][2];
			b2 = Organism[j][3];
			m2 = tan(Organism[j][4]*0.0174532925);
			p2 = Organism[j][5];
			s = Organism[j][6];
			//c = Organism[j];
			scx = pow(10.0,Organism[j][9]); //x-scale factor for data points
			scy = pow(10.0,Organism[j][10]); //y-scale factor for data points
			sc = scy/scx;
			fx = pow(10.0,Organism[j][11]); //Scaling factor for x-errorbars.  fx<=1
			fy = pow(10.0,Organism[j][12]); //  and for y-errorbars
			//Intercept of two lines b1+m1(xint-p1) = b2+m2(xint-p2)
			xint = (b2-m2*p2-b1+m1*p1)/(m1-m2);
			//For BHc2, or RVc2 with theta2>0, need to find where derivative of model=0, to help bracket tangent points  
			if (m1*m2 >= 0) // if slopes both have the same sign, there is no place where dydx=0; use intersection point
							//  Note: this should never happen for BH vs. c2!
			{
				x0 = xint; 
			}
			else // x0 = point where dydx=0
			{
				if(m1<0)
				{
					x0 = (log(-m1)-log(m2)-s*b2+s*m2*p2+s*b1-s*m1*p1)/(s*(m2-m1));
				}
				else
				{
					x0 = (log(-m2)-log(m1)-s*b1+s*m1*p1+s*b2-s*m2*p2)/(s*(m1-m2));
				}
			}
			pt.X[0]=x0;
			y0 = (f)(Organism[j],pt); // y-value of peak of BH vs. c2 curve
			
			for(n=0;n<(signed)SLI.Data[j].size();n++)
			{
				p = &SLI.Data[j][n];
				pt = *p;
				pn = pt;
				xn = p->X[0]; // data point x-value
				sigxn = p->X_PLUS[0]; // data point 1-sigma x-errorbar (assume symmetric Gaussian)
				yn = p->Y; // data point y-value
				//sigyn = p->Y_PLUS; // data point 1-sigma y-errorbar (assume symmetric Gaussian)
				sigxn = sigxn*fx;
				//sigyn = sigyn*fy;
				sigyn = 0.0;
				weight = p->X[1]; // data point weight (should be in second "X column", after data X errorbar
				//cout << "xn = " << xn << " sigxn = " << sigxn << " yn = " << yn << " sigyn = " << sigyn << " weight = " << weight << endl;
				//system("pause");
				ex = sqrt(sigxn*sigxn+slopx*slopx);
				//ey = sqrt(sigyn*sigyn+slopy*slopy);
				ey = 0.0;
				//pn = Data Point with errorbars convolved with slop.
				pn.X_PLUS[0]=ex;
				pn.Y_PLUS=ey; 
				//yxn = (f)(Organism[j],pn); //Value of model curve at x=xn.
				//cout << "xn = " << xn << " ex = " << ex << " yn = " << yn << " ey = " << ey << endl;
				//cout << "x0 = " << x0 << " y0 = " << y0 << " yxn = " << yxn << endl;

				if (yn > y0) // If data point is above peak of curve, there is no y for which y(x)=yn.  PENALIZE!
				{
					chisqr = 1.0E20;
				}

				else // If data point is below peak of curve, there are two points where y(x)=yn
					 // Find the one on the same side of the peak as the data point
				{
					if(xn<x0)
					{
						xt = rtbis(dxy, Organism[j], pn, -3.0, x0, xacc);
					}
					else
					{
						xt = rtbis(dxy, Organism[j], pn, x0, 5.0, xacc);
					}
					chisqr = 2.0*log(ex)+(xn-xt)*(xn-xt)/ex/ex;
				}				
				//if(yn>7)
				//	cout << xn << " " << ex << " " << yn << " " << ey << " " << xt << " " << yt << " " << slope << " " << chisqr << endl;
				chi += weight*chisqr;
			}
		return chi;
		}
	}

#ifdef _WIN32
	_declspec(dllexport) double TRF_BHc2_Y(vector<vector<double> > &Organism)
#else
	double TRF_BHc2_Y(vector<vector<double> > &Organism)
#endif
	{
		//c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s, 
		//c[7]=SlopX, c[8]=Slopy, c[9]=logScaleX, c[10]=logScaleY, c[11]=logfX, c[12]=logfY;

		double slopx, slopy;
		double xn,sigxn,yn,sigyn,weight;
		double ex,ey,e0,arg,fac;
		double chisqr,xt,yt,yxt,slope;
		double xint,yint,b1,m1,p1,b2,m2,p2,s;
		int j,n;
		double chi = 0.0;
		double ax,bx,cx,xmin;
		double f1,f2,x0,y0,x1,x2,yt1,yt2,arg1,arg2,xt1,xt2,slope1,slope2,chisqr1,chisqr2,yxn;
		double e01,e02,fac1,fac2;
		double scx,scy,sc,fx,fy;
		double dx,fmid,f0,xmid,xacc,xdchi2min,dchi2min;
		Point* p;
		Point pt, px1, px2, pn;
		CHIMOD dchi,dchidx,dchidx2,dxy;
		dchi = &DCHI;
		dchidx = &DCHIDX;
		dchidx2 = &DCHIDX2;
		dxy = &DXY;
		pt.X.push_back(0.0);
		xacc = 1.0E-12;
	
		for(j=0;j<(signed)SLI.Data.size();j++)
		{
			slopx = Organism[j][SlopX_Position]; 
			slopy = Organism[j][SlopY_Position];
			//intersection of 2 lines: c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s
			b1 = Organism[j][0];
			m1 = tan(Organism[j][1]*0.0174532925);
			p1 = Organism[j][2];
			b2 = Organism[j][3];
			m2 = tan(Organism[j][4]*0.0174532925);
			p2 = Organism[j][5];
			s = Organism[j][6];
			//c = Organism[j];
			scx = pow(10.0,Organism[j][9]); //x-scale factor for data points
			scy = pow(10.0,Organism[j][10]); //y-scale factor for data points
			sc = scy/scx;
			fx = pow(10.0,Organism[j][11]); //Scaling factor for x-errorbars.  fx<=1
			fy = pow(10.0,Organism[j][12]); //  and for y-errorbars
			//Intercept of two lines b1+m1(xint-p1) = b2+m2(xint-p2)
			//xint = (b2-m2*p2-b1+m1*p1)/(m1-m2);
			//For BHc2, or RVc2 with theta2>0, need to find where derivative of model=0, to help bracket tangent points  
			/*
			if (m1*m2 >= 0) // if slopes both have the same sign, there is no place where dydx=0; use intersection point
							//  Note: this should never happen for BH vs. c2!
			{
				x0 = xint; 
			}
			else // x0 = point where dydx=0
			{
				if(m1<0)
				{
					x0 = (log(-m1)-log(m2)-s*b2+s*m2*p2+s*b1-s*m1*p1)/(s*(m2-m1));
				}
				else
				{
					x0 = (log(-m2)-log(m1)-s*b1+s*m1*p1+s*b2-s*m2*p2)/(s*(m1-m2));
				}
			}
			pt.X[0]=x0;
			y0 = (f)(Organism[j],pt); // y-value of peak of BH vs. c2 curve
			*/
			for(n=0;n<(signed)SLI.Data[j].size();n++)
			{
				p = &SLI.Data[j][n];
				pt = *p;
				//pn = pt;
				xn = p->X[0]; // data point x-value
				//sigxn = p->X_PLUS[0]; // data point 1-sigma x-errorbar (assume symmetric Gaussian)
				yn = p->Y; // data point y-value
				sigyn = p->Y_PLUS; // data point 1-sigma y-errorbar (assume symmetric Gaussian)
				//sigxn = sigxn*fx;
				sigyn = sigyn*fy;
				//sigyn = 0.0;
				sigxn = 0.0;
				weight = p->X[1]; // data point weight (should be in second "X column", after data X errorbar
				//cout << "xn = " << xn << " sigxn = " << sigxn << " yn = " << yn << " sigyn = " << sigyn << " weight = " << weight << endl;
				//system("pause");
				//ex = sqrt(sigxn*sigxn+slopx*slopx);
				ex = 0.0;
				ey = sqrt(sigyn*sigyn+slopy*slopy);
				//ey = 0.0;
				//pn = Data Point with errorbars convolved with slop.
				//pn.X_PLUS[0]=ex;
				//pn.Y_PLUS=ey; 
				//yxn = (f)(Organism[j],pn); //Value of model curve at x=xn.
				//cout << "xn = " << xn << " ex = " << ex << " yn = " << yn << " ey = " << ey << endl;
				//cout << "x0 = " << x0 << " y0 = " << y0 << " yxn = " << yxn << endl;
				xt = xn;
				pt.X[0] = xt;
				yt = (f)(Organism[j],pt);
				chisqr = 2.0*log(sc)+2.0*log(ey)+(yn-yt)*(yn-yt)/ey/ey;
				chi += weight*chisqr;
			}
		return chi;
		}
	}

#ifdef _WIN32
	_declspec(dllexport) double TRF_RVc2(vector<vector<double> > &Organism)
#else
	double TRF_RVc2(vector<vector<double> > &Organism)
#endif
	{
		//c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s, 
		//c[7]=SlopX, c[8]=Slopy, c[9]=logScaleX, c[10]=logScaleY, c[11]=logfX, c[12]=logfY;

		double slopx, slopy;
		double xn,sigxn,yn,sigyn,weight;
		double ex,ey,e0,arg,fac;
		double chisqr,xt,yt,yxt,slope;
		double xint,yint,b1,m1,p1,b2,m2,p2,s;
		int j,n;
		double chi = 0.0;
		double ax,bx,cx,xmin;
		double f1,f2,x0,y0,x1,x2,yt1,yt2,arg1,arg2,xt1,xt2,slope1,slope2,chisqr1,chisqr2,yxn;
		double e01,e02,fac1,fac2;
		double scx,scy,sc,fx,fy;
		double dx,fmid,f0,xmid,xacc,xdchi2min,dchi2min;
		Point* p;
		Point pt, px1, px2, pn;
		CHIMOD dchi,dchidx,dchidx2;
		dchi = &DCHI;
		dchidx = &DCHIDX;
		dchidx2 = &DCHIDX2;
		pt.X.push_back(0.0);
		xacc = 1.0E-6;
	
		/*
		double sumwt = 0.0;
		double sumlnsig = 0.0;
		double sumterm = 0.0;
		double sumlnfac = 0.0;
		double sumfac = 0.0;
		double sumchi = 0.0;
		double sumsig = 0.0;
		double sumr = 0.0;
		double sumr2 = 0.0;
		double sumr4 = 0.0;
		double summ2 = 0.0;
		double sumar = 0.0;
		double sumbr = 0.0;
		double sumax = 0.0;
		double sumbx = 0.0;
		double sumay = 0.0;
		double sumby = 0.0;
		double r,ai,bi;
		*/
		//string filename;
		//char* file;

		//filename = "diagnose.txt";
		//file = const_cast<char*> (filename.c_str());
		//std::ofstream tmp;
		//tmp.open(name);
		//ofstream outfile("diagnose.txt");
		for(j=0;j<(signed)SLI.Data.size();j++)
		{
			slopx = Organism[j][SlopX_Position]; 
			slopy = Organism[j][SlopY_Position];
			//intersection of 2 lines: c[0]=b1, c[1]=t1, c[2]=p1, c[3]=b2, c[4]=t2, c[5]=p2, c[6]=s
			b1 = Organism[j][0];
			m1 = tan(Organism[j][1]*0.0174532925);
			p1 = Organism[j][2];
			b2 = Organism[j][3];
			m2 = tan(Organism[j][4]*0.0174532925);
			p2 = Organism[j][5];
			s = Organism[j][6];
			//c = Organism[j];
			scx = pow(10.0,Organism[j][9]); //x-scale factor for data points
			scy = pow(10.0,Organism[j][10]); //y-scale factor for data points
			sc = scy/scx;
			fx = pow(10.0,Organism[j][11]); //Scaling factor for x-errorbars.  fx<=1
			fy = pow(10.0,Organism[j][12]); //  and for y-errorbars
			//Intercept of two lines b1+m1(xint-p1) = b2+m2(xint-p2)
			xint = (b2-m2*p2-b1+m1*p1)/(m1-m2);
			//For BHc2, or RVc2 with theta2>0, need to find where derivative of model=0, if it exists, to help bracket tangent points  
			if (m1*m2 >= 0) // if slopes both have the same sign, there is no place where dydx=0; set x0 = intersection point of two linear segments
				
			{
				x0 = xint; 
			}
			else // x0 = point where dydx=0
			{
				if(m1<0)
				{
					x0 = (log(-m1)-log(m2)-s*b2+s*m2*p2+s*b1-s*m1*p1)/(s*(m2-m1));
				}
				else
				{
					x0 = (log(-m2)-log(m1)-s*b1+s*m1*p1+s*b2-s*m2*p2)/(s*(m1-m2));
				}
			}
			pt.X[0]=x0;
			y0 = (f)(Organism[j],pt); // y-value of minimum of RV vs. c2 curve (or of intersection point)
			
			for(n=0;n<(signed)SLI.Data[j].size();n++)
			{
				p = &SLI.Data[j][n];
				pt = *p;
				pn = pt;
				xn = p->X[0]; // data point x-value
				sigxn = p->X_PLUS[0]; // data point 1-sigma x-errorbar (assume symmetric Gaussian)
				yn = p->Y; // data point y-value
				sigyn = p->Y_PLUS; // data point 1-sigma y-errorbar (assume symmetric Gaussian)
				sigxn = sigxn*fx;
				sigyn = sigyn*fy;
				weight = p->X[1]; // data point weight (should be in second "X column", after data X errorbar
				//cout << "xn = " << xn << " sigxn = " << sigxn << " yn = " << yn << " sigyn = " << sigyn << " weight = " << weight << endl;
				//system("pause");
				ex = sqrt(sigxn*sigxn+slopx*slopx);
				ey = sqrt(sigyn*sigyn+slopy*slopy);
				//pn = Data Point with errorbars convolved with slop.
				pn.X_PLUS[0]=ex;
				pn.Y_PLUS=ey; 
				yxn = (f)(Organism[j],pn); //Value of model curve at x=xn.
				//cout << "xn = " << xn << " ex = " << ex << " yn = " << yn << " ey = " << ey << endl;
				//cout << "x0 = " << x0 << " y0 = " << y0 << " yxn = " << yxn << endl;

				if (((m1*m2<0.0) && (yn < y0)) || ((m1*m2>=0.0) && (yn < yxn)))
					         // For RVc2, if the curve has a minimum and the data point is below the minimum, 
							 //  Or, if it doesn't have a minimum and the data point is below the curve,
							 //  there is only one tangent point, on the same side of the
							 //  peak as the data point...find it using rtbis on DCHIDX
				{
					//cout << " yn below: " << endl;
					x1 = -3.0;
					x2 = 5.0;
					pt.X[0]=x1;
					f1 = (dchidx)(Organism[j],pt,pn);
					pt.X[0]=x2;
					f2 = (dchidx)(Organism[j],pt,pn);
					//In RVc2 fits with theta2<0, ex>>ey, root of dchidx is often well to the right of c2=5.0.  Expand search range if necessary.
					while(f1*f2>=0.0)
					{
						//x1 = -5.0;
						x2 *= 2.0;
						//pt.X[0]=x1;
						//f1 = (dchidx)(Organism[j],pt,pn);
						pt.X[0]=x2;
						f2 = (dchidx)(Organism[j],pt,pn);
						/*
						if(f1*f2>=0.0)
						{
							cout << "dchidx root not bracketed between -5 and 10!" << endl;
							cout << "xn sigxn yn sigyn ex ey" << endl;
							cout << xn << " " << sigxn << " " << yn << " " << sigyn << " " << ex << " " << ey << endl;
							cout << "b1 " << b1 << endl;
							cout << "t1 " << Organism[j][1] << endl;
							cout << "p1 " << p1 << endl;
							cout << "b2 " << b2 << endl;
							cout << "t2 " << Organism[j][4] << endl;
							cout << "p2 " << p2 << endl;
							cout << "s " << s << endl;
							cout << "SlopX " << slopx << endl;
							cout << "SlopY " << slopy << endl;
							cout << "logsX " << log10(scx) << endl;
							cout << "logsY " << log10(scy) << endl;
							cout << "logfX " << log10(fx) << endl;
							cout << "logfY " << log10(fy) << endl;
							system("pause");
						}
						*/
					}
					if (xn<x0)
					{
						xt = rtbis(dchidx,Organism[j],pn,x1,x2,xacc);
					}
					else
					{
						xt = rtbis(dchidx,Organism[j],pn,x1,x2,xacc);
					}
					pt.X[0]=xt;
					yt = (f)(Organism[j],pt);
					slope = (dfdx)(Organism[j],pt);
					//cout << " xt = " << xt << " yt = " << yt << " slope = " << slope << endl;
				}

				else // If data point is above minimum of curve, or above the curve when it has no minimum, there can be up to 3 tangent points
								   // In either case, 2nd derivative DCHIDX2 will have a single minimum, somewhere near x0.
								   //  Find it using mnbrak and brent
								   // If the minimum is > 0, there is only one tangent point
								   // Otherwise, there may be more than one
				{
					//cout << " yn above:" << endl;
					ax = x0;
					bx = x0+0.1;
					mnbrak(dchidx2,Organism[j],pn,ax,bx,cx);
					xdchi2min = brent(dchidx2,Organism[j],pn,ax,bx,cx);
					pt.X[0]=xdchi2min;
					dchi2min = (dchidx2)(Organism[j],pt,pn);
					//cout << " xdchi2min = " << xdchi2min << " dchi2min = " << dchi2min << endl;
					if (dchi2min >= 0.0) // If dchi2min > 0, there is only one tangent point
						                 //   Find it using rtbis on dchidx
					{
						pt.X[0]=-3.0;
						f1 = (dchidx)(Organism[j],pt,pn);
						pt.X[0]=5.0;
						f2 = (dchidx)(Organism[j],pt,pn);
						if (f1*f2>=0.0)
						{
							cout << "dchidx root not bracketed between -3 and 5!" << endl;
							system("pause");
						}
						if(xn > x0)
						{
							xt = rtbis(dchidx, Organism[j], pn, -3.0, 5.0, xacc);
						}
						else
						{
							xt = rtbis(dchidx, Organism[j], pn, -3.0, 5.0, xacc);
						}
						pt.X[0]=xt;
						yt = (f)(Organism[j],pt);
						slope = (dfdx)(Organism[j],pt);
						//cout << " dchi2min > 0:" << endl;
						//cout << " xt = " << xt << " yt = " << yt << " slope = " << slope << endl;
					}
					else  // If dchi2min < 0, there may be up to three tangent points, one of which is at a local maximimum of dchi
						  // The tangent points, if they exist, are bracketed by zeroes of dchidx2
						  //  Find the bracketing points x1, x2 using rtbis on dchidx2
					{
						//cout << " dchi2min < 0" << endl;
						x1 = rtbis(dchidx2, Organism[j], pn, -3.0, xdchi2min, xacc);
						x2 = rtbis(dchidx2, Organism[j], pn, xdchi2min, 5.0, xacc);
						//cout << " x1 = " << x1 << " x2 = " << x2 << endl;
						// The local maximum of dchi is between x1 and x2...I think.  Ignore it.
						// Find the local minima tangent points < x1 and > x2, if they exist, using rtbis on dchidx
						pt.X[0]=-3.0;
						f1 = (dchidx)(Organism[j],pt,pn);
						pt.X[0]=x1;
						f2 = (dchidx)(Organism[j],pt,pn);
						//cout <<  "f1 = " << f1 << " f2 = " << f2 << endl;
						if (f1*f2 < 0.0)
							xt1 = rtbis(dchidx, Organism[j], pn, -3.0, x1, xacc);
						else
							xt1 = -20.0;
						pt.X[0]=xt1;
						yt1 = (f)(Organism[j],pt);
						slope1 = (dfdx)(Organism[j],pt);
						e01 = sqrt(slope1*slope1*ex*ex + ey*ey);
						arg1 = yn-(yt1+slope1*(xn-xt1));
						if(ex == 0.0)
							fac1 = 1.0;
						else
							fac1 = (slope1*slope1+ey*ey/ex/ex)/sqrt(slope1*slope1+sc*sc*pow(ey/ex,4));
						chisqr1 = arg1*arg1/e01/e01 - 2.0*log(fac1) + 2.0*log(e01);
						
						pt.X[0]=x2;
						f1=(dchidx)(Organism[j],pt,pn);
						pt.X[0]=5.0;
						f2=(dchidx)(Organism[j],pt,pn);
						//cout << " f1 = " << f1 << " f2 = " << f2 << endl;
						if(f1*f2 < 0.0)
							xt2 = rtbis(dchidx, Organism[j], pn, x2, 5.0, xacc);
						else
							xt2 = 20.0;
						pt.X[0]=xt2;
						yt2 = (f)(Organism[j],pt);
						slope2 = (dfdx)(Organism[j],pt);
						e02 = sqrt(slope2*slope2*ex*ex + ey*ey);
						arg2 = yn-(yt2+slope2*(xn-xt2));
						if(ex == 0.0)
							fac2 = 1.0;
						else
							fac2 = (slope2*slope2+ey*ey/ex/ex)/sqrt(slope2*slope2+sc*sc*pow(ey/ex,4));
						chisqr2 = arg2*arg2/e02/e02 - 2.0*log(fac2) + 2.0*log(e02);
						//cout << "xt1 = " << xt1 << " yt1 = " << yt1 << " chisqr1 = " << chisqr1 << endl;
						//cout << "xt2 = " << xt2 << " yt2 = " << yt2 << " chisqr2 = " << chisqr2 << endl;

						//Choose the tangent point that minimizes chisqr
						if(chisqr1<chisqr2)
						{
							xt=xt1;
							yt=yt1;
							slope = slope1;
						}
						else
						{
							xt=xt2;
							yt=yt2;
							slope = slope2;
						}
					}
				}

				e0 = sqrt(slope*slope*ex*ex + ey*ey);
				arg = yn-(yt+slope*(xn-xt));
				if(ex == 0.0)
					fac = 1.0;
				else
					fac = (slope*slope+ey*ey/ex/ex)/sqrt(slope*slope+sc*sc*pow(ey/ex,4));
				chisqr = arg*arg/e0/e0 - 2.0*log(fac) + 2.0*log(e0);
				
				//outfile << xn << " " << ex << " " << yn << " " << ey << " " << xt << " " << yt << " " << weight << " " << slope << " " << chisqr << endl;
				//cout << xn << " " << ex << " " << yn << " " << ey << " " << xt << " " << yt << " " << weight << " " << slope << " " << chisqr << endl;
				chi += weight*chisqr;
				//system("pause");
				/*
				r = ey/ex;
				ai = 4.0*sc*sc*r*r*r/(slope*slope+sc*sc*r*r*r*r)-4.0*r/(slope*slope+r*r);
				bi = 1.0/e0/e0-chisqr/e0/e0;
				sumwt += weight;
				sumlnsig += weight*log(e0);
				sumsig += weight*e0;
				sumr += weight*r;
				sumr2 += weight*r*r;
				sumr4 += weight*r*r*r*r;
				summ2 += weight*slope*slope;
				sumchi += weight*chisqr;
				sumlnfac += weight*log((slope*slope+r*r)/sqrt(slope*slope+sc*sc*r*r*r*r));
				sumfac += weight*(slope*slope+r*r)/sqrt(slope*slope+sc*sc*r*r*r*r);
				sumar += weight*ai;
				sumbr += weight*2.0*slopx*slopy*(1.0-slope*slope/r/r)*bi;
				sumax += -weight*slopy/slopx/slopx*ai;
				sumbx += weight*2.0*slope*slope*slopx*bi;
				sumay += weight/slopx*ai;
				sumby += weight*2.0*slopy*bi;
				*/
			}

		}
		//Return -2*ln(prob). Sign gets changed in BayesianInference.dll
		/*
		cout << "chi = " << chi << endl;
		cout << " N = " << sumwt << endl;
		cout << " SUMlnfac = " << sumlnfac << endl;
		cout << " SUMlnSIG = " << sumlnsig << endl;
		cout << " SUMchi = " << sumchi << endl;
		cout << " <r> = " << sumr/sumwt << endl;
		cout << " <r^2> = " << sumr2/sumwt << endl;
		cout << " <r^4> = " << sumr4/sumwt << endl;
		cout << " <m^2> = " << summ2/sumwt << endl;
		cout << " <SIG> = " << sumsig/sumwt << endl;
		cout << " <fac> = " << sumfac/sumwt << endl;
		cout << " <Ar> = " << sumar/sumwt << endl;
		cout << " <Br> = " << sumbr/sumwt << endl;
		cout << " <Ax> = " << sumax/sumwt << endl;
		cout << " <Bx> = " << sumbx/sumwt << endl;
		cout << " <Ay> = " << sumay/sumwt << endl;
		cout << " <By> = " << sumby/sumwt << endl;
		system("pause");
		*/
		//outfile.close();
		//system("pause");
		return chi;
	}
}