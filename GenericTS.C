#include<iostream>
using namespace std;

//function for calculating integral of histogram from xmin to xmax
double h_integral(TH1F *h,double xmin,double xmax){

TAxis *axis = h->GetXaxis();
  int bmin = axis->FindBin(xmin); //in your case xmin=-1.5
  int bmax = axis->FindBin(xmax); //in your case xmax=0.8
  double integral = h->Integral(bmin,bmax);
  integral -= (h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin)))/axis->GetBinWidth(bmin);
  integral -= (h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax))/axis->GetBinWidth(bmax);

return integral;

}


//test statistic function defined as the negative log of likelihood ratios
double qmu(double mu,double s,double b, int n, int nobs){

double Ldata_qmu= ((pow(mu*s+b,n))*exp(-(mu*s+b)));

double Ldata_qmuhat = ((pow(nobs,n))*exp(-(nobs)));

double TS_mu = -2*log(Ldata_qmu/Ldata_qmuhat);

return TS_mu;
}



void GenericTS(){

TRandom3 r;

//Data Card details
///////////////////

double mu[7]={1.95,1.35,1.05,1.20,1.31,1.3438,1.3436};
int nobs=45;//number of observed events

double s=20.; //signal rate
double b=30.; //bkg rate


int nsb=0;//random variable for toy generation assuming signal+bkg 
int nb=0; //random variable for toy generation assuming bkg 
int Ntoy=s+b;//s+b rate
///////////////////
cout<<"Make global fit of asimov data"<<endl;

for(int mui=0;mui<7;mui++){//looping over different signal strength


//test statistic histogram distribution from toy s+b only data
TH1F *t1 =new TH1F("t1","t1",330,-30,40);
t1->SetLineColor(2);
t1->SetLineWidth(2);

//test statistic histogram distribution from toy b only data
TH1F *t2 =new TH1F("t2","t2",330,-30,40);
t2->SetLineColor(4);
t2->SetLineWidth(2);




for(int i=0;i<100000;i++){
nsb=r.Poisson(Ntoy);
nb=r.Poisson(b);

t1->Fill(qmu(1,s,b,nsb,nobs));
t2->Fill(qmu(1,s,b,nb,nobs));





}


char Legname1[100];
TLegend *leg_1D[5];
for(int k0=0;k0<5;k0++){
sprintf(Legname1,"leg_1D%i",k0);
leg_1D[k0]=new TLegend(0.2,0.6,0.30,0.8);
leg_1D[k0]->SetTextFont(62);
leg_1D[k0]->SetLineColor(1);
leg_1D[k0]->SetLineStyle(1);
leg_1D[k0]->SetLineWidth(3);
leg_1D[k0]->SetFillColor(0);
leg_1D[k0]->SetFillStyle(1001);
leg_1D[k0]->SetShadowColor(0);
leg_1D[k0]->SetDrawOption(0);
leg_1D[k0]->SetBorderSize(0);
leg_1D[k0]->SetTextSize(0.03);
}

leg_1D[0]->AddEntry(t1,"f(q_{#mu=1}|#mu = 1) , using signal+bkg toy","l");
leg_1D[0]->AddEntry(t2,"f(q_{#mu=1}|#mu = 0) , using bkg toy","l");


t1->Scale(1/t1->Integral());

t2->Scale(1/t2->Integral());
t1->GetYaxis()->SetRangeUser(0,1.5*max(t1->GetMaximum(),t2->GetMaximum()));
t1->Draw();
t2->Draw("SAME");

leg_1D[0]->Draw();

double qmu_obs=qmu(mu[mui],s,b,nobs,nobs);//observed value of test statistic

double p_mu=h_integral(t1,qmu_obs,40);// CLsb
double One_minus_pb=h_integral(t2,qmu_obs,40); //CLb


double CLs =p_mu/One_minus_pb;


cout<<"At r = "<<mu[mui]<<"   :   "<<"q_mu = "<<qmu_obs<<"       q_A = "<<" ---"<<"      CLsb = "<<p_mu<<"     CLb = "<<One_minus_pb<<"          CLs :  "<<CLs<<endl;
if(mui !=6){
delete t1;
delete t2;
}
}//looping over different signal strength


}
