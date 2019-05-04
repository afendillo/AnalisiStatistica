#define Pi TMath::Pi() 
#define e TMath::E()
#define NRand 1e7
#define N 30
#define bin 100
#define xk pow(2 , 1)
using namespace std;
//-----------------------------------------------------------------------------------
//--------------------- Delete old root objects root --------------------------------
//-----------------------------------------------------------------------------------
void Reset(){
	//Recupera la lista di tutti i canvas e li cancella
	TSeqCollection* canvases = gROOT->GetListOfCanvases();
	TIter next(canvases);
	while(TCanvas *c = (TCanvas*)next())
	{
		delete c;
	}
	//Cancella tutti gli oggetti creati
	//Necessario per evitare "warning <TROOT::Append>: Replacing existing TH1: <Name> (Potential memory leak)"
	gROOT->DeleteAll();
	return;
}
//----------------------------------------------------------------------
//----------------------Funzione Fattoriale -----------------------------
//----------------------------------------------------------------------
long double Fattoriale(long double n)
{
  if (n == 0)
    return 1;
  else
    return(n * Fattoriale(n-1));
}
//-----------------------------------------------------------------------
//--------------------Nepero: Definizione Limite ------------------------
//-----------------------------------------------------------------------
double NeperoLimite(long double n)
{
	return pow((n+1)/n ,n);
}
//-----------------------------------------------------------------------
//-----------------------Sviluppo in serie ------------------------------
//-----------------------------------------------------------------------
double NeperoTaylor( int It , int n=0)
{
	if(n>It) return 0;
	return (1/Fattoriale(n))+NeperoTaylor(It,n+1);
}
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
float NeperoTaylorFloat( int It , int n=0)
{
	if(n>It) return 0;
	return (1/(float)Fattoriale(n))+NeperoTaylor(It,n+1);
}
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// double NeperoTaylorVar( int It , int n=0)
// {
	// if(n>It) return 0;
	// return (pow(xk,n)/Fattoriale(n))+NeperoTaylorVar(It,n+1);
// }
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// double NeperoTaylorSim(int It)
// {
	// double sum=0;
	// for(int i =0 ; i<=It ; ++i) sum += 1.0/((double)Fattoriale((double) i));
	// return sum;
// }
//-----------------------------------------------------------------------
//----------------Sviluppo in serie con Min.Comun.Den--------------------
//-----------------------------------------------------------------------
double NeperoTaylorV2(int n)
{
	double Somma=1 , Var=1;
	for(int i=0 ; i<n ; i++)
	{
		for(int x=n ; x>i ; x--) Var=Var*x;
		Somma+=Var;
		Var=1;
	}
	return (Somma/Fattoriale(n));
}
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
float NeperoTaylorV2Float(int n)
{
	float Somma=1, Var=1, esp=0;
	for(int i=0 ; i<n ; i++)
	{
		for(int x=n ; x>i ; x--) Var=Var*x;
		Somma+=Var;
		Var=1;
	}
	return (Somma/Fattoriale(n));
}
//-----------------------------------------------------------------------
//--------- Algoritmo di Horner: raccoglimento di fattori----------------
//-----------------------------------------------------------------------
float HornerFloat(int It, float n=0)
{
	if(It==0) return 1;
	if(n==It) return 0;
	if(n<1) return 2+HornerFloat(It,n+1);
	return (1./(n+1))*(1+HornerFloat(It,n+1));
}
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
double Horner(int It, double n=0)
{
	if(It==0) return 1;
	if(n==It) return 0;
	if(n<1) return 2+Horner(It,n+1);
	return (1./(n+1))*(1+Horner(It,n+1));
}
//-----------------------------------------------------------------------
//------------------- Sviluppo in frazioni continue ---------------------
//-----------------------------------------------------------------------
double Fcontinua(double c , int n=N)
{
	static int counter=0;
	if(counter<1) 
	{
		//cout<<counter<<" ";
		counter++;
		return (int)c+Fcontinua(c , n);
	}

	if(counter>n)
	{
		counter=0;
		return 0;
	}
	double a=1/((c-(int)c));
	//cout<<setprecision(8)<<c*1.<<" "<<1.*int(c)<<"\n";
	//cout<<counter<<" ";
	counter++;
	return 1/((int)a+Fcontinua(a , n));
}

float FcontinuaFloat(float c , int n=N)
{
	static int counter=0;
	if(counter<1) 
	{
		//cout<<counter<<" ";
		counter++;
		return (int)c+FcontinuaFloat(c , n);
	}

	if(counter>n)
	{
		counter=0;
		return 0;
	}
	float a=1/((c-(int)c));
	//cout<<setprecision(8)<<c*1.<<" "<<1.*int(c)<<"\n";
	//cout<<counter<<" ";
	counter++;
	return 1/((int)a+FcontinuaFloat(a , n));
}
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// void PartialFracton(int n , double a)
// {
	//vector <double> FunzContinua;
	// double c=a;
	// for (int i=0 ; i<n ; i++)
	// {
		//FunzContinua.push_back((int)c);
		// cout<<(int)c<<" ";
		// if((c-(int)c)==0) break;
		// c=1/((c-(int)c));
		
	// }
	// cout<<"\n";

	//return 1/((int)c+PartialFracton(n , c));
	
// }
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// void ScientificNotation(double n)
// {	
	// int Iterator=0;
	// double Ten=0;
	// string Sign;
	
	// if(n<1) {Ten=10; Sign="-";}
	// else {Ten=0.1; Sign="+";}
	
	// for(int i=0; (int)n<1 || (int)n>10 ; i++)
	// {
		// n=n*Ten;
		// Iterator=i;
	// }
	// cout<<n<<"E"<<Sign;
	// if (Iterator<10) cout<<"0";
	// cout<<Iterator;
// }
//-----------------------------------------------------------------------
//----------------------------- Main ------------------------------------
//-----------------------------------------------------------------------
void Nepero()

{
	
	Reset();
	//cout<<setprecision(40)
	cout<<"\n";
	cout<<"e = "<<e<<"\n";
	cout<<"\n---------------------------- Limite ------------------------------------------------------------\n\n";
	cout<<"Nepero Float Limit : "<<(float)NeperoLimite(N)<<" Variazione : "<<abs((float)NeperoLimite(N)-e)<<"\n";
	cout<<"Nepero Double Limit : "<<NeperoLimite(N)<<" Variazione : "<<abs(NeperoLimite(N)-e)<<"\n";
	
	cout<<"\n--------------------------- Sviluppo di Taylor --------------------------------------------------\n\n";
	cout<<"Nepero Float Taylor ordine "<<N <<" : " <<NeperoTaylorFloat(N)<<" Variazione : "<<1e6*abs(NeperoTaylorFloat(N)-e)<<"E-06\n";
	cout<<"Nepero Double Taylor ordine "<<N <<" : "<<NeperoTaylor(N)<<" Variazione : "<<1e6*(e-NeperoTaylor(N))<<"E-06\n";
	
	cout<<"\n------------------------------Sviluppo di Taylor McD ----------------------------------------------\n\n";
	cout<<"Nepero Float TaylorV2 ordine "<<N <<" : "<<NeperoTaylorV2Float(N)<<" Variazione : "<<1e6*abs(NeperoTaylorV2Float(N)-e)<<"E-06\n";
	cout<<"Nepero Double TaylorV2 ordine "<<N <<" : "<<NeperoTaylorV2(N)<<" Variazione : "<<(NeperoTaylorV2(N)-e)<<"\n";
	
	cout<<"\n------------------------------Sviluppo di Taylor Horner ----------------------------------------------"<<"\n\n";
	cout<<"Nepero Horner float ordine "<<N <<" : "<<HornerFloat(N)<<" Variazione : "<<1e6*abs(HornerFloat(N)-e)<<"E-06"<<"\n";
	cout<<"Nepero Horner double ordine "<<N <<" : "<<Horner(N)<<" Variazione : "<<1e6*abs(Horner(N)-e)<<"E-06"<<"\n";
	
	cout<<"\n------------------------------Sviluppo in frazioni continue ----------------------------------------------\n";
	cout<<"Nepero Frazioni continue ordine "<<N <<" : "<<FcontinuaFloat(e,N)<<" Variazione : "<<1e6*abs(FcontinuaFloat(e,N)-e)<<"E-06\n\n";
	cout<<"Nepero Frazioni continue ordine "<<N <<" : "<<Fcontinua(e,N)<<" Variazione : "<<1e6*abs(Fcontinua(e,N)-e)<<"E-06\n";
	cout<<"\n";
	
	int line=9;
	double width=1;
	
	auto mg = new TMultiGraph();
	auto gr1 = new TGraph(); gr1->SetLineColor(kBlue); gr1->SetLineStyle(line); gr1->SetLineWidth(width);
	auto gr2 = new TGraph(); gr2->SetLineColor(kRed); gr2->SetLineStyle(line); gr2->SetLineWidth(width);
	auto gr3 = new TGraph(); gr3->SetLineColor(kGreen); gr3->SetLineStyle(line); gr3->SetLineWidth(width);
	auto gr4 = new TGraph(); gr4->SetLineColor(kOrange+2); gr4->SetLineStyle(line); gr4->SetLineWidth(width);
	auto gr5 = new TGraph(); gr5->SetLineColor(kViolet); gr5->SetLineStyle(line); gr5->SetLineWidth(width);
	auto gr6 = new TGraph(); gr6->SetLineColor(kBlack); gr6->SetLineStyle(line); gr6->SetLineWidth(width);
	auto gr7 = new TGraph(); gr7->SetLineColor(kAzure); gr7->SetLineStyle(line); gr7->SetLineWidth(width);
	auto gr8 = new TGraph(); gr8->SetLineColor(kMagenta); gr8->SetLineStyle(line); gr8->SetLineWidth(width);
	
	for (int i=0 ; i<N ; i++)
	{
		//if(abs(NeperoTaylor(i)-e)<abs(NeperoTaylorSim(i)-e)) 
		//{	
		///	cout<<i<<" ";
		//}
		gr1->SetPoint(i,i,e-(float)NeperoLimite(i));
		gr2->SetPoint(i,i,abs(e-NeperoLimite(i)));
		gr3->SetPoint(i,i,abs(e-NeperoTaylorFloat(i)));

		if (abs(e-NeperoTaylor(i))>0) gr4->SetPoint(i,i,abs(e-NeperoTaylor(i)));
		else gr4->SetPoint(i,i,(double)(1e-15));

		gr5->SetPoint(i,i,abs(e-NeperoTaylorV2Float(i)));

		if (abs(e-NeperoTaylorV2(i))>0) gr6->SetPoint(i,i,(double)(1e-15)+abs(e-NeperoTaylorV2(i)));
		else gr6->SetPoint(i,i,(double)(1e-15));

		if (abs(e-FcontinuaFloat(e,i))>0) gr7->SetPoint(i,i,abs(e-FcontinuaFloat(e,i)));
		else gr7->SetPoint(i,i,(double)(1e-15));

		if (abs(e-Fcontinua(e,i))>0) gr8->SetPoint(i,i,(double)(1e-15)+abs(e-Fcontinua(e,i)));
		else gr8->SetPoint(i,i,(double)(1e-15));		
		
		//else cout<<"Simo"<<"\n";
		//cout<<e<<"-"<<NeperoTaylor(i)<<"="<<(e-NeperoTaylor(i))<<"\n";
	}
	cout<<"\n";
	
	mg->Add(gr1); gr1->SetTitle("LimFloat"); gr1->SetMarkerStyle(26); gr1->SetMarkerColor(kBlue);
	mg->Add(gr2); gr2->SetTitle("LimDouble"); gr2->SetMarkerStyle(kPlus); gr2->SetMarkerColor(kRed);
	mg->Add(gr3); gr3->SetTitle("TaylorFloat"); gr3->SetMarkerStyle(kCircle); gr3->SetMarkerColor(kGreen);
	mg->Add(gr4); gr4->SetTitle("TaylorDouble"); gr4->SetMarkerStyle(kStar); gr4->SetMarkerColor(kOrange+2);
	mg->Add(gr5); gr5->SetTitle("TaylorV2Float"); gr5->SetMarkerStyle(kMultiply); gr5->SetMarkerColor(kViolet);
	mg->Add(gr6); gr6->SetTitle("TaylorV2Double");gr6->SetMarkerStyle(kOpenSquare); gr6->SetMarkerColor(kBlack);
	mg->Add(gr7); gr7->SetTitle("FContineFloat"); gr7->SetMarkerStyle(kStar); gr7->SetMarkerColor(kAzure);
	mg->Add(gr8); gr8->SetTitle("FContinueDouble");gr8->SetMarkerStyle(kPlus); gr8->SetMarkerColor(kMagenta);
	
	mg->SetTitle("Algoritmo Nepero; Iterazioni; |e-e_{alg}|");
	
	//mg->GetHistogram()->GetYaxis()->SetRangeUser(1e-16,e-1);
	
	double x=0.65 , y=0.2;
	
	auto legend = new TLegend(x-0.05,y,x+0.25 ,y+0.2);
	legend->SetTextSize(0.024);
	legend->SetTextFont(42);
//	legend->SetHeader("Algoritmi","C"); // option "C" allows to center the header
	legend->AddEntry(gr1,"Nepero Limite Float","lp");
	legend->AddEntry(gr2,"Nepero Limite Double","lp");
	
	legend->AddEntry(gr3,"Nepero Taylor Float","lp");
	legend->AddEntry(gr5,"Nepero TaylorV2 Float","lp");
	legend->AddEntry(gr7,"Nepero Fraz.Continue Float","lp");
	
	legend->AddEntry(gr4,"Nepero Taylor Double","lp");
	legend->AddEntry(gr6,"Nepero TaylorV2 Double","lp");
	legend->AddEntry(gr8,"Nepero Fraz.Continue Double","lp");	
	
	TCanvas* c1 = new TCanvas();
	c1->SetGrid();
	c1 -> SetLogy();
	mg->Draw("apc");
	//mg->SetTitleSize(50);
	mg->GetXaxis()->SetTitleSize(0.055);mg->GetXaxis()->SetTitleOffset(0.7);
  mg->GetYaxis()->SetTitleSize(0.055);mg->GetYaxis()->SetTitleOffset(0.85);
	legend->Draw();
	// gr2->Draw("pc same");
	// gr3->Draw("pc same");
	// gr4->Draw("pc same");
	// gr5->Draw("pc same");
	// gr6->Draw("pc same");

	// for (int i=0 ; i<gr7->GetN(); i++)
	// {
	// 	cout<<gr8->GetY()[i]<<endl;
	// }
	
	return;
}