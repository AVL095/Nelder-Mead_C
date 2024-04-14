#include <string.h>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <vector> 
#include <map> 
#include <numeric>
#include <functional>
#include <ctime>


int Nelder_Mead(int nk,vector<double>&x0,double step,int lim,double tol,double(*func)(vector<double>)){
    double alpha,gamma,rho,sigma,fr,fe,fc,max_diff,sum,q;
    int i,j,counter,kk;

     //init_simplex

    alpha = 1.0; gamma = 2.0; rho = 0.5; sigma = 0.5;
    
    vector<vector<double>>simplex(nk+1,vector<double>(nk,0.0));
    vector<double>direction(nk,0.0);
    vector<double>f_values(nk+1);
    vector<int>order(nk+1);
    vector<double>centroid(nk,0.0);
    vector<double>xr(nk);
    vector<double>xe(nk);
    vector<double>xc(nk);

    simplex[0]=x0;

    for (i=0;i<nk;++i) {
        direction[i]=step;
        simplex[i+1]=x0;
        for (j=0;j<nk;++j) {
            simplex[i+1][j]+=direction[j];
        }
    }
    
//#########################Начало итераций##########################################
    counter=0;
 
    for (kk=0;kk<lim;kk++) {
        counter++;
        for(i=0;i<nk+1;++i) f_values[i]=func(simplex[i]); //evaluate(оценка симплекса)
        iota(order.begin(),order.end(),0);
        sort(order.begin(),order.end(),[&](int i, int j){return f_values[i]<f_values[j];});

        vector<vector<double>>new_simplex(nk+1,vector<double>(nk));
        vector<double>new_f_values(nk+1);

        for (i=0;i<nk+1;++i) {
            new_simplex[i]=simplex[order[i]];new_f_values[i]=f_values[order[i]];
        }
        simplex=move(new_simplex);f_values=move(new_f_values);

    //#################### Условие выхода по достижению точности###################
        //max_diff = 0.0;
        //for (double fi:f_values) max_diff=max(max_diff,abs(fi-f_values[0]));
        //if (max_diff<=tol) break;
          q=func(simplex[0]);
          if(q<=tol) break;
    //#############################################################################    
 
    //Центроид всех точек, кроме худшей

    for (i=0;i<nk;++i)
    centroid[i]=accumulate(simplex.begin(),simplex.end()-1,0.0,[i](double sum,const vector<double>&point){return sum+point[i];})/nk;

        for(i=0;i<nk;++i) xr[i]=centroid[i]+alpha*(centroid[i]-simplex.back()[i]); //reflect (отражение)
        fr=func(xr);
//##########################################################################################
    if (f_values[0] <= fr && fr < f_values[nk - 1]) {
            simplex.back()=xr;
        } else if (fr<f_values[0]) {
           for(i=0;i<nk;++i) xe[i]=centroid[i]+gamma*(xr[i]-centroid[i]); //expand(растяжение)
            fe=func(xe);
            simplex.back()=fe<fr ? xe:xr;
        } else {
            for(i=0;i<nk;++i) xc[i]=centroid[i]+rho*(simplex.back()[i]-centroid[i]); //contract(сжатие)

            fc=func(xc);
           //###############  
            if (fc<f_values.back()) {
                simplex.back() = xc;
             } else {
                //shrink(уменьшение размера)
                   for (i=1;i<nk+1;++i) {
                      for(j=0;j<nk;++j) simplex[i][j]=simplex[0][j]+sigma*(simplex[i][j]-simplex[0][j]);
                   }
            }
          //#################
    }
//#############################################################################################
 } //end kk

//########################
    x0=simplex[0];
    return counter;
}

//######################################################

double  NormalMinFunction(vector<double>xsimpl) {
    double s1,s2,s3,s4,z,psi,p,d,c1,c2;
    int i,kx;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;
    if (xsimpl[0]<=0) return 10000;
    if (xsimpl[1]<=0) return 10000;

    for (i=0;i<nesm.n;i++) {
            z = (nesm.x[i]-xsimpl[0])/xsimpl[1];
            d = spi * exp(-z * z / 2.);
            p = norm_cdf(z);
            psi = d / (1. - p);
            s1 +=(1.-nesm.r[i])*(nesm.x[i]-xsimpl[0]);
            s2 += (1.-nesm.r[i])*pow(nesm.x[i]-xsimpl[0],2);
            s3 += nesm.r[i]*psi;
            s4 += nesm.r[i]*psi*z;
            kx+=1-nesm.r[i];
    }
    c1=s1+xsimpl[1]*s3;
    c2=s2+pow(xsimpl[1],2)*(s4-kx);
    z=c1*c1+c2*c2;
    return z;
}


//#################MLE Weibull Minimized Function#######################

double WeibullMinFunction(vector<double>xsimpl) {
 double s1,s2,s3,z,b,c;
 int i,k;
 if (xsimpl[0]<=0) return(10000000.);
   s1=0;s2=0;s3=0;k=0;
   b=xsimpl[0];
 for(i=0;i<nesm.n;i++) {
   k+=(1-nesm.r[i]);
   s1+=pow(nesm.x[i],b);
 }
   c=s1/k;

for(i=0;i<nesm.n;i++) {
  z=(pow(nesm.x[i],b))/c;
  s3+=z*log(z);
  s2+=(1-nesm.r[i])*log(z);
}
 c=s3-s2-k;
 return c*c;
}


//###############################################################

void MLE_Normal(string ff) {
    
    int i,j,k,kx,icount,kp,nx,lim;
    string s1;
    double **v,*fcum,*ycum,*xplow,*xpup,*zp,*p,tlow,tup,*xp;
    double cp,cko,q,eps,beta,delta,step,*t,tp,z;
//#################################################################
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");
    inp >> s1; //text
    inp >>nesm.n; //sample size
    inp >> s1;
    inp >> beta;
    inp >> s1;
    inp >> step;
    inp >> s1;
    inp >> eps;
    inp >> s1;
    inp >> lim;
    inp >> s1; //text
    for (i = 0; i <nesm.n; i++) {inp>>z;nesm.x.push_back(z);}
    inp >> s1; //text
    for (i=0;i<nesm.n; i++) {inp>>j;nesm.r.push_back(j);}
    inp >> s1;
    inp>>kp;
    p=new double[kp];
    zp=new double[kp];
    xp=new double[kp];
    inp >> s1;
    for (i = 0; i < kp; i++)  inp >> p[i];
    inp.close();
//###########################################################
     nx=2;   
     vector<double>xsimpl;
     v = new double* [nx];
     double *xq;
     int *rq;
     xq=new double [nesm.n];
     rq=new int [nesm.n];

     for(i=0;i<nesm.n;i++) {xq[i]=nesm.x[i];rq[i]=nesm.r[i];}

    for (i = 0; i < nx; i++) v[i] = new double[nx];
    for (i = 0; i < nx; i++) {
        for (j = 0; j < nx; j++) v[i][j] = 0;
    }

    cp=0;cko=0;k=0;
    for(i=0;i<nesm.n;i++) {
        k+=(1-nesm.r[i]); // количество наблюдений
        cp+=(1-nesm.r[i])*nesm.x[i];
        cko+=(1-nesm.r[i])*nesm.x[i]*nesm.x[i];
    }
    cp/=k; //выборочное среднее по наблюдениям
    cko=sqrt((cko-cp*cp*k)/(k-1)); //выборочное ско по наблюдениям
    q=0; 
    xsimpl.push_back(cp);xsimpl.push_back(cko);
    v[0][0]=1.;v[1][1]=0.5;v[0][1]=0.;v[1][0]=0.;

     q=0;
     //eps=1.0e-15;lim=500;
   
    if(k!=nesm.n) {
     icount=Nelder_Mead(nx,xsimpl,step,lim,eps,NormalMinFunction);
     q=NormalMinFunction(xsimpl);
     CovMatrixMleN(nesm.n,xq,rq, xsimpl[0],xsimpl[1],v);
    }

    kx=kp;
    int nc;
    nc=4;
    int m[]={k,kx,kx,kx};

    t=new double[kx];
    xplow=new double[kx];
    xpup=new double[kx];
     
    for (i = 0; i < kx; i++) {
      zp[i]=norm_ppf(p[i]);
      xp[i]=xsimpl[0]+zp[i]*xsimpl[1];
      delta = zp[i] * sqrt(nesm.n);
      if(k==nesm.n) {
          t[i]=nct_ppf(beta,nesm.n-1,delta);
      }
      else {
           lmtaprn(beta,nesm.n-1,v[0][0],v[1][1],v[0][1],delta,tlow,tup);
           xpup[i]=xsimpl[0]+tup*xsimpl[1]/sqrt(nesm.n); 
           xplow[i]=xsimpl[0]+tlow*xsimpl[1]/sqrt(nesm.n); 
      }
     }
         if(k==nesm.n) {
           for (i = 0; i < kx; i++) xpup[i]=xsimpl[0]+t[i]*xsimpl[1]/sqrt(nesm.n); 
           for (i = 0; i < kx; i++) xplow[i]=xsimpl[0]-t[kx-i-1]*xsimpl[1]/sqrt(nesm.n);
         }      

//########################################################
    vector<double>xx;vector<double>y; 

    fcum = new double[nesm.n];ycum = new double[nesm.n];
    cum(nesm.n,xq,rq,k,fcum,ycum);

    for (i = 0; i < k; i++) {
       xx.push_back(ycum[i]);
       y.push_back( 5.0 + norm_ppf(fcum[i]));
    }

    for (i = 0; i < kx; i++) {
       xx.push_back(xplow[i]);
       y.push_back(zp[i]+5.0);
    }

    for (i = 0; i < kx; i++) {
       xx.push_back(xp[i]);
       y.push_back(zp[i]+5.0);
    }
    
    for (i = 0; i < kx; i++) {
       xx.push_back(xpup[i]);
       y.push_back(zp[i]+5.0);
    }
 
    ofstream out1("Out/" + ff + ".xout");
    out1<<nc<<endl;
    for(i=0;i<nc;i++) out1<<m[i]<<" ";
    out1<<endl;

    int km=0;
    for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<xx[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    }
     
     km=0;
     for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<y[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    } 
     out1.close();

 //############################################################### 

 out << "Method:" << ff << "\n";
 out << "n=" <<nesm.n << "\n";
 out << "X" << "\n";
 for (i = 0; i <nesm.n; i++)  out <<nesm.x[i] << " , ";
 out << "\n";
 out << "R" << "\n";
 for (i = 0; i<nesm.n; i++)  out <<nesm.r[i] << " , ";
 out << "\n";
 out << "cp*=" <<setprecision(12) << fixed <<cp << "\n";
 out << "cko*="<<setprecision(12) << fixed <<cko << "\n";
 out << "Q="<<setprecision(12) << q << "\n";
 out << "icount="<<icount<<endl;
 out << "cp="<<setprecision(12) << fixed << xsimpl[0]<<"\n";
 out << "cko="<<setprecision(12) << fixed <<xsimpl[1]<<"\n";
 out << "P" << "\n";
 for (i = 0; i < kx; i++) out << p[i] << " ; ";
 out << "\n";
 out << "Xplow" << "\n";
 for (i = 0; i < kx; i++) out << xplow[i] << " ; ";
 out << "\n";
 out << "Xp" << "\n";
 for (i = 0; i < kx; i++) out << xp[i] << " ; ";
 out << "\n";
 out << "Xpup" << "\n";
 for (i = 0; i < kx; i++) out << xpup[i] << " ; ";
 out << "\n";
 out << "v11=" << v[0][0] << "\n";
 out << "v12=" << v[0][1] << "\n";
 out << "v21=" << v[1][0] << "\n";
 out << "v22=" << v[1][1] << "\n";
 out.close();

 delete[] fcum, ycum, v,xplow,xpup,zp,p,xp;
 nesm.r.clear();nesm.x.clear();xsimpl.clear();xx.clear();y.clear();

}

//#########################################################################

void MLE_Weibull(string ff) {
 int i,j,k,kx,kp,nx,icount,lim;
 string s1;
 double **v, *fcum, *ycum,*xplow,*xpup,*p,*zp,*xp,tlow,tup;
 double cp, cko,eps, q,s,b,c,aw,sw,z,step,beta;

 vector<double>logx;
 vector<double>xsimpl;
//##################################################################		   
 ifstream inp("Inp/" + ff + ".inp");
 inp >> s1;
 inp >>nesm.n;
 inp >> s1;
 inp >> beta;
 inp >> s1;
 inp >> step;
 inp >> s1;
 inp >> eps;
 inp >> s1;
 inp >> lim;
 inp >> s1;
 for(i=0;i<nesm.n;i++) {inp>>z;nesm.x.push_back(z);}
 inp >> s1;
 for(i= 0;i<nesm.n;i++) {inp>>j;nesm.r.push_back(j);}
 inp >> s1;
 inp>> kp;
 p=new double[kp];
 zp=new double[kp];
 xp=new double[kp];
 inp >> s1;
 for (i = 0; i < kp; i++)  inp >> p[i];
 inp.close();
//###############################################################
 cp = 0;cko=0;k=0;
 for(i=0;i<nesm.n;i++) logx.push_back(log(nesm.x[i]));
 for(i=0;i<nesm.n;i++) {
   k +=(1-nesm.r[i]); // количество наблюдений
   cp+=(1-nesm.r[i])*logx[i];
   cko+=(1-nesm.r[i])*logx[i]*logx[i];
   }
 cp/=k; //выборочное среднее по наблюдениям
 cko=sqrt((cko-cp*cp*k)/(k-1)); //выборочное ско по наблюдениям
 v = new double* [2];
 for (i = 0; i < 2; i++) v[i] = new double[2];
 for (i = 0; i < 2; i++) {
   for (j = 0; j < 2; j++) v[i][j] = 0;
 }

   double *xq;
   int *rq;
   xq=new double [nesm.n];
   rq=new int [nesm.n];

   for(i=0;i<nesm.n;i++) {xq[i]=nesm.x[i];rq[i]=nesm.r[i];}

 nx=1;q=0;xsimpl.push_back(1./cko);

 icount=Nelder_Mead(nx,xsimpl,step,lim,eps,WeibullMinFunction);
 q=WeibullMinFunction(xsimpl);

 b=xsimpl[0];
 s = 0.;
 for(i=0;i<nesm.n;i++) s+=pow(nesm.x[i],b);
 c=pow(s/k,1/b);
 CovMatrixMleW(nesm.n,xq,rq,c,b,v);
 kx=kp;
 sw=1./b;aw=log(c);

 int nc;
 nc = 4;
 int m[]={k,kx,kx,kx};
 vector<double>xx;
 vector<double>y; 
 xplow=new double[kx];
 xpup=new double[kx];
 for (i = 0; i < kx; i++) {
   zp[i] = log(log(1 / (1 - p[i])));
   xp[i]=aw+sw*zp[i];
   nctWeibull(beta,nesm.n-1,v[0][0],v[1][1],v[0][1],zp[i]*sqrt(nesm.n),tlow,tup);
   xplow[i]=aw+tlow*sw/sqrt(nesm.n);
   xpup[i]=aw+tup*sw/sqrt(nesm.n);
}

//########################################################

    fcum = new double[nesm.n];
    ycum = new double[nesm.n];
    cum(nesm.n,xq,rq,k,fcum,ycum);
   
 for (i = 0; i < k; i++) {
        xx.push_back(log(ycum[i]));
        y.push_back(log(log(1./(1 - fcum[i])))+5.);
    }

 for (i = 0; i < kx; i++) {
       xx.push_back(xplow[i]);
       y.push_back(zp[i]+5.);
 }

 for (i = 0; i < kx; i++) {
       xx.push_back(xp[i]);
       y.push_back(zp[i]+5.);
 }
                    
  for (i = 0; i < kx; i++) {
      xx.push_back(xpup[i]);
      y.push_back(zp[i]+5.);
  }

    ofstream out1("Out/" + ff + ".xout");
    out1<<nc<<endl;
    for(i=0;i<nc;i++) out1<<m[i]<<" ";
    out1<<endl;

    int km=0;
    for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<xx[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    }
     
     km=0;
     for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<y[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    } 
     out1.close();

//###########################################################################
  ofstream out("Out/" + ff + ".out");
  out << "Method:" << ff << "\n";
  out << "n=" <<nesm.n << "\n";
  out << "X" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<nesm.x[i] << " , ";
  out << "\n";
  out << "log(X)" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<logx[i] << " , ";
  out << "\n";
  out << "R" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<nesm.r[i] << " , ";
  out << "\n";
  out << "cp*="<<setprecision(12)<< fixed<<cp << "\n";
  out << "cko*="<<setprecision(12)<< fixed<<cko << "\n";
  out << "Q=" << q << "\n";
  out <<"icount="<<icount<<endl;
  out << "b="<<setprecision(12)<< fixed <<b<<"\n";
  out << "c="<<setprecision(12) << fixed <<c<< "\n";
  out << "aw="<<setprecision(12)<< fixed <<aw<<"\n";
  out << "sw="<<setprecision(12) << fixed <<sw<< "\n";
  out << "P" << "\n";
  for (i = 0; i < kx; i++) out << p[i] << " ; ";
  out << "\n";
  out << "Xplow" << "\n";
  for (i = 0; i < kx; i++) out << xplow[i] << " ; ";
  out << "\n";
  out << "Xp" << "\n";
  for (i = 0; i < kx; i++) out << xp[i] << " ; ";
  out << "\n";
  out << "Xpup" << "\n";
  for (i = 0; i < kx; i++) out << xpup[i] << " ; ";
  out << "\n";
  out << "v11=" << v[0][0] << "\n";
  out << "v12=" << v[0][1] << "\n";
  out << "v21=" << v[1][0] << "\n";
  out << "v22=" << v[1][1] << "\n";
  out.close();
 
  logx.clear();nesm.r.clear();xsimpl.clear();nesm.x.clear();xx.clear();y.clear();
  delete[] fcum,ycum,v,xplow,xpup,p,zp,xp;
  }




int Nelder_Mead(int nk,vector<double>&x0,double step,int lim,double tol,double(*func)(vector<double>)){
    double alpha,gamma,rho,sigma,fr,fe,fc,max_diff,sum,q;
    int i,j,counter,kk;

     //init_simplex

    alpha = 1.0; gamma = 2.0; rho = 0.5; sigma = 0.5;
    
    vector<vector<double>>simplex(nk+1,vector<double>(nk,0.0));
    vector<double>direction(nk,0.0);
    vector<double>f_values(nk+1);
    vector<int>order(nk+1);
    vector<double>centroid(nk,0.0);
    vector<double>xr(nk);
    vector<double>xe(nk);
    vector<double>xc(nk);

    simplex[0]=x0;

    for (i=0;i<nk;++i) {
        direction[i]=step;
        simplex[i+1]=x0;
        for (j=0;j<nk;++j) {
            simplex[i+1][j]+=direction[j];
        }
    }
    
//#########################Начало итераций##########################################
    counter=0;
 
    for (kk=0;kk<lim;kk++) {
        counter++;
        for(i=0;i<nk+1;++i) f_values[i]=func(simplex[i]); //evaluate(оценка симплекса)
        iota(order.begin(),order.end(),0);
        sort(order.begin(),order.end(),[&](int i, int j){return f_values[i]<f_values[j];});

        vector<vector<double>>new_simplex(nk+1,vector<double>(nk));
        vector<double>new_f_values(nk+1);

        for (i=0;i<nk+1;++i) {
            new_simplex[i]=simplex[order[i]];new_f_values[i]=f_values[order[i]];
        }
        simplex=move(new_simplex);f_values=move(new_f_values);

    //#################### Условие выхода по достижению точности###################
        //max_diff = 0.0;
        //for (double fi:f_values) max_diff=max(max_diff,abs(fi-f_values[0]));
        //if (max_diff<=tol) break;
          q=func(simplex[0]);
          if(q<=tol) break;
    //#############################################################################    
 
    //Центроид всех точек, кроме худшей

    for (i=0;i<nk;++i)
    centroid[i]=accumulate(simplex.begin(),simplex.end()-1,0.0,[i](double sum,const vector<double>&point){return sum+point[i];})/nk;
/*
     for (i=0;i<nk;i++) {
        sum= 0.0;
        for (j=0;j<nk;j++) sum+=simplex[j][i];
        centroid[i]=sum/nk;
    }
*/
        for(i=0;i<nk;++i) xr[i]=centroid[i]+alpha*(centroid[i]-simplex.back()[i]); //reflect (отражение)
        fr=func(xr);
//##########################################################################################
    if (f_values[0] <= fr && fr < f_values[nk - 1]) {
            simplex.back()=xr;
        } else if (fr<f_values[0]) {
           for(i=0;i<nk;++i) xe[i]=centroid[i]+gamma*(xr[i]-centroid[i]); //expand(растяжение)
            fe=func(xe);
            simplex.back()=fe<fr ? xe:xr;
        } else {
            for(i=0;i<nk;++i) xc[i]=centroid[i]+rho*(simplex.back()[i]-centroid[i]); //contract(сжатие)

            fc=func(xc);
           //###############  
            if (fc<f_values.back()) {
                simplex.back() = xc;
             } else {
                //shrink(уменьшение размера)
                   for (i=1;i<nk+1;++i) {
                      for(j=0;j<nk;++j) simplex[i][j]=simplex[0][j]+sigma*(simplex[i][j]-simplex[0][j]);
                   }
            }
          //#################
    }
//#############################################################################################
 } //end kk

//########################
    x0=simplex[0];
    return counter;
}

//######################################################

double  NormalMinFunction(vector<double>xsimpl) {
    double s1,s2,s3,s4,z,psi,p,d,c1,c2;
    int i,kx;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;
    if (xsimpl[0]<=0) return 10000;
    if (xsimpl[1]<=0) return 10000;

    for (i=0;i<nesm.n;i++) {
            z = (nesm.x[i]-xsimpl[0])/xsimpl[1];
            d = spi * exp(-z * z / 2.);
            p = norm_cdf(z);
            psi = d / (1. - p);
            s1 +=(1.-nesm.r[i])*(nesm.x[i]-xsimpl[0]);
            s2 += (1.-nesm.r[i])*pow(nesm.x[i]-xsimpl[0],2);
            s3 += nesm.r[i]*psi;
            s4 += nesm.r[i]*psi*z;
            kx+=1-nesm.r[i];
    }
    c1=s1+xsimpl[1]*s3;
    c2=s2+pow(xsimpl[1],2)*(s4-kx);
    z=c1*c1+c2*c2;
    return z;
}


//#################MLE Weibull Minimized Function#######################

double WeibullMinFunction(vector<double>xsimpl) {
 double s1,s2,s3,z,b,c;
 int i,k;
 if (xsimpl[0]<=0) return(10000000.);
   s1=0;s2=0;s3=0;k=0;
   b=xsimpl[0];
 for(i=0;i<nesm.n;i++) {
   k+=(1-nesm.r[i]);
   s1+=pow(nesm.x[i],b);
 }
   c=s1/k;

for(i=0;i<nesm.n;i++) {
  z=(pow(nesm.x[i],b))/c;
  s3+=z*log(z);
  s2+=(1-nesm.r[i])*log(z);
}
 c=s3-s2-k;
 return c*c;
}


//###############################################################

void MLE_Normal(string ff) {
    
    int i,j,k,kx,icount,kp,nx,lim;
    string s1;
    double **v,*fcum,*ycum,*xplow,*xpup,*zp,*p,tlow,tup,*xp;
    double cp,cko,q,eps,beta,delta,step,*t,tp,z;
//#################################################################
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");
    inp >> s1; //text
    inp >>nesm.n; //sample size
    inp >> s1;
    inp >> beta;
    inp >> s1;
    inp >> step;
    inp >> s1;
    inp >> eps;
    inp >> s1;
    inp >> lim;
    inp >> s1; //text
    for (i = 0; i <nesm.n; i++) {inp>>z;nesm.x.push_back(z);}
    inp >> s1; //text
    for (i=0;i<nesm.n; i++) {inp>>j;nesm.r.push_back(j);}
    inp >> s1;
    inp>>kp;
    p=new double[kp];
    zp=new double[kp];
    xp=new double[kp];
    inp >> s1;
    for (i = 0; i < kp; i++)  inp >> p[i];
    inp.close();
//###########################################################
     nx=2;   
     vector<double>xsimpl;
     v = new double* [nx];
     double *xq;
     int *rq;
     xq=new double [nesm.n];
     rq=new int [nesm.n];

     for(i=0;i<nesm.n;i++) {xq[i]=nesm.x[i];rq[i]=nesm.r[i];}

    for (i = 0; i < nx; i++) v[i] = new double[nx];
    for (i = 0; i < nx; i++) {
        for (j = 0; j < nx; j++) v[i][j] = 0;
    }

    cp=0;cko=0;k=0;
    for(i=0;i<nesm.n;i++) {
        k+=(1-nesm.r[i]); // количество наблюдений
        cp+=(1-nesm.r[i])*nesm.x[i];
        cko+=(1-nesm.r[i])*nesm.x[i]*nesm.x[i];
    }
    cp/=k; //выборочное среднее по наблюдениям
    cko=sqrt((cko-cp*cp*k)/(k-1)); //выборочное ско по наблюдениям
    q=0; 
    xsimpl.push_back(cp);xsimpl.push_back(cko);
    v[0][0]=1.;v[1][1]=0.5;v[0][1]=0.;v[1][0]=0.;

     q=0;
     //eps=1.0e-15;lim=500;
   
    if(k!=nesm.n) {
     icount=Nelder_Mead(nx,xsimpl,step,lim,eps,NormalMinFunction);
     q=NormalMinFunction(xsimpl);
     CovMatrixMleN(nesm.n,xq,rq, xsimpl[0],xsimpl[1],v);
    }

    kx=kp;
    int nc;
    nc=4;
    int m[]={k,kx,kx,kx};

    t=new double[kx];
    xplow=new double[kx];
    xpup=new double[kx];
     
    for (i = 0; i < kx; i++) {
      zp[i]=norm_ppf(p[i]);
      xp[i]=xsimpl[0]+zp[i]*xsimpl[1];
      delta = zp[i] * sqrt(nesm.n);
      if(k==nesm.n) {
          t[i]=nct_ppf(beta,nesm.n-1,delta);
      }
      else {
           lmtaprn(beta,nesm.n-1,v[0][0],v[1][1],v[0][1],delta,tlow,tup);
           xpup[i]=xsimpl[0]+tup*xsimpl[1]/sqrt(nesm.n); 
           xplow[i]=xsimpl[0]+tlow*xsimpl[1]/sqrt(nesm.n); 
      }
     }
         if(k==nesm.n) {
           for (i = 0; i < kx; i++) xpup[i]=xsimpl[0]+t[i]*xsimpl[1]/sqrt(nesm.n); 
           for (i = 0; i < kx; i++) xplow[i]=xsimpl[0]-t[kx-i-1]*xsimpl[1]/sqrt(nesm.n);
         }      

//########################################################
    vector<double>xx;vector<double>y; 

    fcum = new double[nesm.n];ycum = new double[nesm.n];
    cum(nesm.n,xq,rq,k,fcum,ycum);

    for (i = 0; i < k; i++) {
       xx.push_back(ycum[i]);
       y.push_back( 5.0 + norm_ppf(fcum[i]));
    }

    for (i = 0; i < kx; i++) {
       xx.push_back(xplow[i]);
       y.push_back(zp[i]+5.0);
    }

    for (i = 0; i < kx; i++) {
       xx.push_back(xp[i]);
       y.push_back(zp[i]+5.0);
    }
    
    for (i = 0; i < kx; i++) {
       xx.push_back(xpup[i]);
       y.push_back(zp[i]+5.0);
    }
 
    ofstream out1("Out/" + ff + ".xout");
    out1<<nc<<endl;
    for(i=0;i<nc;i++) out1<<m[i]<<" ";
    out1<<endl;

    int km=0;
    for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<xx[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    }
     
     km=0;
     for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<y[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    } 
     out1.close();

 //############################################################### 

 out << "Method:" << ff << "\n";
 out << "n=" <<nesm.n << "\n";
 out << "X" << "\n";
 for (i = 0; i <nesm.n; i++)  out <<nesm.x[i] << " , ";
 out << "\n";
 out << "R" << "\n";
 for (i = 0; i<nesm.n; i++)  out <<nesm.r[i] << " , ";
 out << "\n";
 out << "cp*=" <<setprecision(12) << fixed <<cp << "\n";
 out << "cko*="<<setprecision(12) << fixed <<cko << "\n";
 out << "Q="<<setprecision(12) << q << "\n";
 out << "icount="<<icount<<endl;
 out << "cp="<<setprecision(12) << fixed << xsimpl[0]<<"\n";
 out << "cko="<<setprecision(12) << fixed <<xsimpl[1]<<"\n";
 out << "P" << "\n";
 for (i = 0; i < kx; i++) out << p[i] << " ; ";
 out << "\n";
 out << "Xplow" << "\n";
 for (i = 0; i < kx; i++) out << xplow[i] << " ; ";
 out << "\n";
 out << "Xp" << "\n";
 for (i = 0; i < kx; i++) out << xp[i] << " ; ";
 out << "\n";
 out << "Xpup" << "\n";
 for (i = 0; i < kx; i++) out << xpup[i] << " ; ";
 out << "\n";
 out << "v11=" << v[0][0] << "\n";
 out << "v12=" << v[0][1] << "\n";
 out << "v21=" << v[1][0] << "\n";
 out << "v22=" << v[1][1] << "\n";
 out.close();

 delete[] fcum, ycum, v,xplow,xpup,zp,p,xp;
 nesm.r.clear();nesm.x.clear();xsimpl.clear();xx.clear();y.clear();

}

//#########################################################################

void MLE_Weibull(string ff) {
 int i,j,k,kx,kp,nx,icount,lim;
 string s1;
 double **v, *fcum, *ycum,*xplow,*xpup,*p,*zp,*xp,tlow,tup;
 double cp, cko,eps, q,s,b,c,aw,sw,z,step,beta;

 vector<double>logx;
 vector<double>xsimpl;
//##################################################################		   
 ifstream inp("Inp/" + ff + ".inp");
 inp >> s1;
 inp >>nesm.n;
 inp >> s1;
 inp >> beta;
 inp >> s1;
 inp >> step;
 inp >> s1;
 inp >> eps;
 inp >> s1;
 inp >> lim;
 inp >> s1;
 for(i=0;i<nesm.n;i++) {inp>>z;nesm.x.push_back(z);}
 inp >> s1;
 for(i= 0;i<nesm.n;i++) {inp>>j;nesm.r.push_back(j);}
 inp >> s1;
 inp>> kp;
 p=new double[kp];
 zp=new double[kp];
 xp=new double[kp];
 inp >> s1;
 for (i = 0; i < kp; i++)  inp >> p[i];
 inp.close();
//###############################################################
 cp = 0;cko=0;k=0;
 for(i=0;i<nesm.n;i++) logx.push_back(log(nesm.x[i]));
 for(i=0;i<nesm.n;i++) {
   k +=(1-nesm.r[i]); // количество наблюдений
   cp+=(1-nesm.r[i])*logx[i];
   cko+=(1-nesm.r[i])*logx[i]*logx[i];
   }
 cp/=k; //выборочное среднее по наблюдениям
 cko=sqrt((cko-cp*cp*k)/(k-1)); //выборочное ско по наблюдениям
 v = new double* [2];
 for (i = 0; i < 2; i++) v[i] = new double[2];
 for (i = 0; i < 2; i++) {
   for (j = 0; j < 2; j++) v[i][j] = 0;
 }

   double *xq;
   int *rq;
   xq=new double [nesm.n];
   rq=new int [nesm.n];

   for(i=0;i<nesm.n;i++) {xq[i]=nesm.x[i];rq[i]=nesm.r[i];}

 nx=1;q=0;xsimpl.push_back(1./cko);

 icount=Nelder_Mead(nx,xsimpl,step,lim,eps,WeibullMinFunction);
 q=WeibullMinFunction(xsimpl);

 b=xsimpl[0];
 s = 0.;
 for(i=0;i<nesm.n;i++) s+=pow(nesm.x[i],b);
 c=pow(s/k,1/b);
 CovMatrixMleW(nesm.n,xq,rq,c,b,v);
 kx=kp;
 sw=1./b;aw=log(c);

 int nc;
 nc = 4;
 int m[]={k,kx,kx,kx};
 vector<double>xx;
 vector<double>y; 
 xplow=new double[kx];
 xpup=new double[kx];
 for (i = 0; i < kx; i++) {
   zp[i] = log(log(1 / (1 - p[i])));
   xp[i]=aw+sw*zp[i];
   nctWeibull(beta,nesm.n-1,v[0][0],v[1][1],v[0][1],zp[i]*sqrt(nesm.n),tlow,tup);
   xplow[i]=aw+tlow*sw/sqrt(nesm.n);
   xpup[i]=aw+tup*sw/sqrt(nesm.n);
}

//########################################################

    fcum = new double[nesm.n];
    ycum = new double[nesm.n];
    cum(nesm.n,xq,rq,k,fcum,ycum);
   
 for (i = 0; i < k; i++) {
        xx.push_back(log(ycum[i]));
        y.push_back(log(log(1./(1 - fcum[i])))+5.);
    }

 for (i = 0; i < kx; i++) {
       xx.push_back(xplow[i]);
       y.push_back(zp[i]+5.);
 }

 for (i = 0; i < kx; i++) {
       xx.push_back(xp[i]);
       y.push_back(zp[i]+5.);
 }
                    
  for (i = 0; i < kx; i++) {
      xx.push_back(xpup[i]);
      y.push_back(zp[i]+5.);
  }

    ofstream out1("Out/" + ff + ".xout");
    out1<<nc<<endl;
    for(i=0;i<nc;i++) out1<<m[i]<<" ";
    out1<<endl;

    int km=0;
    for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<xx[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    }
     
     km=0;
     for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<y[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    } 
     out1.close();

//###########################################################################
  ofstream out("Out/" + ff + ".out");
  out << "Method:" << ff << "\n";
  out << "n=" <<nesm.n << "\n";
  out << "X" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<nesm.x[i] << " , ";
  out << "\n";
  out << "log(X)" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<logx[i] << " , ";
  out << "\n";
  out << "R" << "\n";
  for (i = 0; i <nesm.n; i++)  out <<nesm.r[i] << " , ";
  out << "\n";
  out << "cp*="<<setprecision(12)<< fixed<<cp << "\n";
  out << "cko*="<<setprecision(12)<< fixed<<cko << "\n";
  out << "Q=" << q << "\n";
  out <<"icount="<<icount<<endl;
  out << "b="<<setprecision(12)<< fixed <<b<<"\n";
  out << "c="<<setprecision(12) << fixed <<c<< "\n";
  out << "aw="<<setprecision(12)<< fixed <<aw<<"\n";
  out << "sw="<<setprecision(12) << fixed <<sw<< "\n";
  out << "P" << "\n";
  for (i = 0; i < kx; i++) out << p[i] << " ; ";
  out << "\n";
  out << "Xplow" << "\n";
  for (i = 0; i < kx; i++) out << xplow[i] << " ; ";
  out << "\n";
  out << "Xp" << "\n";
  for (i = 0; i < kx; i++) out << xp[i] << " ; ";
  out << "\n";
  out << "Xpup" << "\n";
  for (i = 0; i < kx; i++) out << xpup[i] << " ; ";
  out << "\n";
  out << "v11=" << v[0][0] << "\n";
  out << "v12=" << v[0][1] << "\n";
  out << "v21=" << v[1][0] << "\n";
  out << "v22=" << v[1][1] << "\n";
  out.close();
 
  logx.clear();nesm.r.clear();xsimpl.clear();nesm.x.clear();xx.clear();y.clear();
  delete[] fcum,ycum,v,xplow,xpup,p,zp,xp;
  }


//###########################################################################

void MLE_UpDown(string ff) {

    int nx,i,j,kx,nn,kp,icount,kfail,knon,ksigne,lim,km;
    double **v,*xp,*xplow,*xpup,zp,tplow,tpup;
    double q,eps,stepx,cp,s,bint,delta,d,s1,s2,s3,s4,z,beta,step;
    string st;
   
    vector<double>xsimpl;
    vector<double>px;
    vector<int>nfailure;
    vector<int>nnon;

//##################################################

    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");
    inp>>st;
    inp>>nesm.klevel;
    nx=2;
    inp >> st;
    inp >> beta;
    inp >> st;
    inp >> step;
    inp >> st;
    inp >> eps;
    inp >> st;
    inp >> lim;
    inp>>st;
    for (i=0;i<nesm.klevel;i++) {inp>>z;nesm.x.push_back(z);}
    inp>>st;
    for (i=0;i<nesm.klevel;i++) {inp>>j;nesm.nsample.push_back(j);}
    inp>>st;
    for (i=0;i<nesm.klevel;i++) {inp>>j;nfailure.push_back(j);}
    inp>>st;
    inp>>kp;
    inp>>st;
    for (i=0;i<kp;i++) {inp>>z;px.push_back(z);}
    inp.close();

//###################################################


    if(ff=="UpDown_Normal") nesm.ts="Normal";
    if(ff=="UpDown_Weibull") {
        for (i=0;i<nesm.klevel;i++) nesm.x[i]=log(nesm.x[i]);
        nesm.ts="Weibull";
    }

    //Приближенный расчет оценок среднего и среднего квадратичного отклонения (cp,s)
    //Dixon W. Т., Mood A. М. J. Amer. Statist. Ass., v. 43, 1948, p. 109.

   
    kfail=0;knon=0;ksigne=1;
    for (i=0;i<nesm.klevel;i++) {
        kfail+=nfailure[i];
        j=nesm.nsample[i]-nfailure[i];
        nnon.push_back(j);
        knon+=j;
    }
    if (kfail<knon) ksigne=-1;  //Расчет ведут по разрушенным образцам
    s1=0;s2=0;s3=0;
    for (i=0;i<nesm.klevel;i++) {
        if (ksigne==-1) {
            s1+=i*nfailure[i];
            s2+=i*i*nfailure[i];
            s3+=nfailure[i];
        }
        if (ksigne==1) {
            s1+=i*nnon[i];
            s2+=i*i*nnon[i];
            s3+=nnon[i];
        }
    }
    d=nesm.x[1]-nesm.x[0];
    cp=nesm.x[0]+d*(s1/s3+ksigne*0.5); //Оценка среднего
    s4=(s3*s2-s1*s1)/(s3*s3);
    s=1.62*d*(s4+0.029); //Оценка ско

    nn=0;
    for (i=0;i<nesm.klevel;i++) {
        z=double(nfailure[i])/double(nesm.nsample[i]);
        nesm.p.push_back(z);
        nn+=nesm.nsample[i];
    }

    //Расчет методом максимального правдоподобия
    sm.cko=s;
    //eps=1.0e-15;lim=500;
    q=0;
    xsimpl.push_back(cp);xsimpl.push_back(s);

    icount=Nelder_Mead(nx,xsimpl,step,lim,eps,UpDownMinimize);
    q=UpDownMinimize(xsimpl);

    for(i=0;i<nesm.klevel;i++) {
        if(nesm.x[i]>xsimpl[0]) {
            stepx=(nesm.x[i]-nesm.x[i-1])/xsimpl[1];
            bint=(xsimpl[0]-nesm.x[i-1])/xsimpl[1];
            break;
        }
    }

    v=new double*[nx];
    for (i = 0; i < nx; i++) v[i] = new double[nx];
    for (i = 0; i < nx; i++) {
        for (j = 0; j < nx; j++) v[i][j] = 0;
    }

    CovMatrixUpDown(nesm.ts,stepx,bint,v);

    int nc;
    nc=4;
    int m []={nesm.klevel,kp,kp,kp};
  
    xp = new double[kp];
    xplow = new double[kp];
    xpup = new double[kp];
    vector<double>xx;
    vector<double>y; 

    kx = kp;
    for (i = 0; i < kp; i++) {
        if (nesm.ts == "Normal") zp=norm_ppf(px[i]);
        if (nesm.ts == "Weibull") zp= log(log(1 / (1 - px[i])));
        xp[i] = xsimpl[0]+zp*xsimpl[1];
        delta=zp*sqrt(nn);
        lmtaprn(beta,nn,v[0][0],v[1][1],v[0][1],delta,tplow,tpup);
        xpup[i] = xsimpl[0] + xsimpl[1] * tpup/sqrt(nn);
        xplow[i] = xsimpl[0] + xsimpl[1] * tplow/sqrt(nn);
    }


//######################################################

    for (i = 0; i < nesm.klevel; i++) {
        xx.push_back(nesm.x[i]);
        y.push_back(nesm.p[i]);
    }

     for (i = 0; i < kp; i++) {
       xx.push_back(xplow[i]);
       y.push_back(px[i]);
     }

     for (i = 0; i < kp; i++) {
      xx.push_back(xp[i]);
      y.push_back(px[i]);
     }
     for (i = 0; i < kp; i++) {
       xx.push_back(xpup[i]);
       y.push_back(px[i]);
     }

//########################################################

    ofstream out1("Out/"+ff+".xout");
    out1<<nc<<endl;
    for(i=0;i<nc;i++) out1<<m[i]<<" ";
    out1<<endl;

    km=0;
    for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<xx[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    }
     
     km=0;
     for(i=0;i<nc;i++) {
      for(j=0;j<m[i];j++) out1<<y[j+km]<<" ";
      km+=m[i];
      out1<<endl;
    } 
     out1.close();

 //############################################################### 

    out << "Method:" << ff << "\n";
    out << "k="<<nesm.klevel<<"\n";
    out << "X" << "\n";
    for (i = 0; i < nesm.klevel; i++)  out <<nesm.x[i] << " , ";
    out << "\n";
    out << "n" << "\n";
    for (i = 0; i <nesm.klevel;i++)  out <<nesm.nsample[i]<<" , ";
    out << "\n";
    out << "nfailure" << "\n";
    for (i=0;i<nesm.klevel;i++) out<<nfailure[i]<<" , ";
    out << "\n";
    out << "An approximate calculation is possible if>0.3=" << s4 << endl;
    out << "cp*="<< fixed << setprecision(12) << cp << "\n";
    out << "cko*="<< fixed << setprecision(12) << s << "\n";
    out << "Q="<< fixed << setprecision(12) << q << "\n";
    out << "icount=" << icount << endl;
    out << "x[0]=" << fixed << setprecision(12) << xsimpl[0] << "\n";
    out << "x[1]=" << fixed << setprecision(12) << xsimpl[1] << "\n";
    out << "P" << "\n";
    for (i = 0; i < kp; i++) out << px[i] << " ; ";
    out << "\n";
    out << "Xplow" << "\n";
    for (i = 0; i < kp; i++) out << xplow[i] << " ; ";
    out << "\n";
    out << "Xp" << "\n";
    for (i = 0; i < kp; i++) out << xp[i] << " ; ";
    out << "\n";
    out << "Xpup" << "\n";
    for (i = 0; i < kp; i++) out << xpup[i] << " ; ";
    out << "\n";
    out << "v11=" << v[0][0] << "\n";
    out << "v12=" << v[0][1] << "\n";
    out << "v22=" << v[1][1] << "\n";


    out.close();
    nesm.x.clear();nesm.p.clear();nesm.nsample.clear();px.clear();
    xx.clear();y.clear();nfailure.clear();nnon.clear();
    delete[] v,xp,xpup,xplow;
}


//##############MLE Up-down Minimize Function Normal,Weibull#################################################################

 double UpDownMinimize(vector<double>xsimpl) {

  double z,p,d,s1,s2,fiz;
  int i;

  if((xsimpl[1]<=0) || (xsimpl[1]>1.5*sm.cko)) return(10000000.);

  s1=0;s2=0;
  for(i=0;i<nesm.klevel;i++) {
    z=(nesm.x[i]-xsimpl[0])/xsimpl[1];
    if(nesm.ts=="Normal") {
      p=norm_cdf(z);
      d=norm_pdf(z);
   }
    if(nesm.ts=="Weibull") {
     p=1.-exp(-exp(z));
     d=exp(z-(exp(z)));
    }
  if(p<=0 || p>=1) return(10000000.0);
   fiz=(nesm.p[i]-p)*nesm.nsample[i]*d/(p*(1.-p));
   s1+=fiz;
   s2+=fiz*z;
 }
 return(s1*s1+s2*s2);
}

