#define Pi TMath::Pi() 
#define Massa 0.776
#define Gamma 0.149

#define NSig 2000
#define NBack 5000
#define bin 85
#define A 8
#define B 5
#define C 3
#define xmax 1.5
#define xmin 0
#define Range xmax-xmin
#define PI 3.1415

using namespace std;


//Inverse function of ERF
float myErfInv(float x){
	float tt1, tt2, lnx, sgn;
	sgn = (x < 0) ? -1.0f : 1.0f;

	x = (1 - x)*(1 + x);        // x = 1 - x*x;
	lnx = logf(x);
	tt1 = 2/(Pi*0.147) + 0.5f * lnx;
	tt2 = 1/(0.147) * lnx;

	return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

//Breit-Wigner Function
double BreitWigner ( double *x , double *p) {
	
	double y = 1/PI*(Range)*p[2]/bin*(p[0]/2)/(pow(x[0]-p[1],2)+p[0]*p[0]/4);
	return y; //p[0] = gamma; p[1] = massa; p[2] = Nsig
}

//Breit-Wigner FixFunction
double BreitWignerFix ( double *x , double *p) {
	
	double y = 1/PI*(Range)*NSig/bin*(p[0]/2)/(pow(x[0]-p[1],2)+p[0]*p[0]/4);
	return y; //p[0] = gamma; p[1] = massa; NSig is fixed
}

//Poly function
double pol2 (double *x, double *par ) {

	double y = par[2]*x[0]*x[0] + par[1]*x[0] + par[0];
	return y;
}

 // Sum of background and signal function
double fitFunction(double *x, double *p) {
      return pol2(x,p) + BreitWigner(x,&p[3]);
   }

 // Sum of background and signal function
double fitFunctionfix(double *x, double *p) {
      return pol2(x,p) + BreitWignerFix(x,&p[3]);
   }
//Constant PDFunction a-b interval
double cFunc(double x) {
		return ROOT::Math::uniform_pdf  (x,0,1,0);
							
}	

void FitMethods(char stamp='Q')
{	
	//Recupera la lista di tutti i canvas e li cancella
	int count=0;
	TSeqCollection* canvases = gROOT->GetListOfCanvases();
	TIter next(canvases);
	while(TCanvas *c = (TCanvas*)next())
	{
		delete c;
		count++;
	}
	if(stamp!='Q') cout<<count<<" Canvases Deleted"<<endl;
	//Cancella tutti gli oggetti creati
	//Necessario per evitare "warning <TROOT::Append>: Replacing existing TH1: <Name> (Potential memory leak)"
	gROOT->DeleteAll();
	
	//Informazioni statistiche da stampare
	gStyle->SetOptFit(1111);

//----------------------------------------------------------------------------------------------------------------------	
//--------------------------------- Inverse Function Method for Breit-Wigner -------------------------------------------	
//----------------------------------------------------------------------------------------------------------------------	


	double Max=Pi/2 , Min=-Pi/2 , random=0;
	TRandom3 *RandGenerator = new TRandom3(time(0));
	TRandom3 *Gauss = new TRandom3((unsigned)time( NULL ));
	TF1 *breit = new TF1("breit",BreitWigner, Massa-Gamma*Range, Massa+Gamma*Range , 1);
	breit->SetParameter(0, 1);
	
	TH1D *BreitHisto = new TH1D ("BreitHisto" , "Least Squares" , bin ,xmin, xmax);
	BreitHisto->GetXaxis()->SetTitle("Massa [GeV]");
	BreitHisto->GetYaxis()->SetTitle("Conteggi");
	BreitHisto->SetFillColorAlpha(kGreen, 0.30);

	TH1D *BreitHisto1 = new TH1D ("BreitHisto1" , "Least Squares" , bin ,xmin, xmax);
	BreitHisto1->GetXaxis()->SetTitle("Massa [GeV]");
	BreitHisto1->GetYaxis()->SetTitle("Conteggi");
	BreitHisto1->SetFillColorAlpha(kGreen, 0.30);

	TH1D *BreitHisto2 = new TH1D ("BreitHisto2" , "Modified Least Squares" , bin ,xmin, xmax);
	BreitHisto2->GetXaxis()->SetTitle("Massa [GeV]");
	BreitHisto2->GetYaxis()->SetTitle("Conteggi");
	BreitHisto2->SetFillColorAlpha(kGreen, 0.30);

	TH1D *BreitHisto3 = new TH1D ("BreitHisto3" , "Modified Least Squares" , bin ,xmin, xmax);
	BreitHisto3->GetXaxis()->SetTitle("Massa [GeV]");
	BreitHisto3->GetYaxis()->SetTitle("Conteggi");
	BreitHisto3->SetFillColorAlpha(kGreen, 0.30);

	TH1D *BreitHisto4 = new TH1D ("BreitHisto4" , "Maximum Likelihood" , bin ,xmin, xmax);
	BreitHisto4->GetXaxis()->SetTitle("Massa [GeV]");
	BreitHisto4->GetYaxis()->SetTitle("Conteggi");
	BreitHisto4->SetFillColorAlpha(kGreen, 0.30);
	
	for (int i = 0 ; i<NSig; i++)
	{
		random=Massa+Gamma/2*tan(RandGenerator->Rndm()*(Max-Min)+Min-TMath::ATan(2*Massa/Gamma));
		random = random + Gauss->Gaus(0,0.01);
		BreitHisto->Fill(random);
		BreitHisto1->Fill(random);
		BreitHisto2->Fill(random);
		BreitHisto3->Fill(random);
		BreitHisto4->Fill(random);	
	}
//----------------------------------------------------------------------------------------------------------------------	
//--------------------------------- Hit or Miss Poly Background --------------------------------------------------------	
//----------------------------------------------------------------------------------------------------------------------		
	double frandx , randx , randy;
// First we need to evaluate the maximum of the poly. we have 3 options:
	double ymax1 = A*xmax*xmax + B*xmax + C;
	double ymax2 = A*xmin*xmin + B*xmin + C;
	double ymax3 = A*((-B/(2*A))*(-B/(2*A))) + B*(-B/(2*A)) + C;

	double ymax = max(ymax1,ymax2);

	ymax = max(ymax,ymax3);
	
	for (int i = 0; i<NBack; i++)
	{ 
// Then we generate a couple of random variables uniformly

		randx = RandGenerator->Rndm()*(xmax-xmin) + xmin;
		randy = RandGenerator->Rndm()*ymax;
		frandx = A*randx*randx + B*randx + C;

//If randy is less than or equals f(randx) we accept the background event, else we try again
	
		if (randy <= frandx)
		{			
			randx = randx + Gauss->Gaus(0,0.01);
			BreitHisto->Fill(randx);
			BreitHisto1->Fill(randx);
			BreitHisto2->Fill(randx);
			BreitHisto3->Fill(randx);
			BreitHisto4->Fill(randx);
		}
		else
			i=i-1;
	}
	

	TF1  *fitfunc = new TF1("fitfunc",fitFunction,0 , 1.5 , 6);
	fitfunc->SetParameter(0, C);
	fitfunc->SetParName(0, "C");
	fitfunc->SetParameter(1, B);
	fitfunc->SetParName(1, "B");
	fitfunc->SetParameter(2, A);
	fitfunc->SetParName(2, "A");
	fitfunc->SetParameter(3, Gamma);
	fitfunc->SetParName(3, "Gamma");
	fitfunc->SetParameter(4, Massa);
	fitfunc->SetParName(4, "Massa");
	fitfunc->SetParameter(5, NSig);
	fitfunc->SetParName(5, "Segnale");

	TF1  *fitfuncfix = new TF1("fitfuncfix",fitFunctionfix,0 , 1.5 , 5);
	fitfuncfix->SetParameter(0, C);
	fitfuncfix->SetParName(0, "C");
	fitfuncfix->SetParameter(1, B);
	fitfuncfix->SetParName(1, "B");
	fitfuncfix->SetParameter(2, A);
	fitfuncfix->SetParName(2, "A");
	fitfuncfix->SetParameter(3, Gamma);
	fitfuncfix->SetParName(3, "Gamma");
	//fitfuncfix->FixParameter(3, Gamma);
	fitfuncfix->SetParameter(4, Massa);
	fitfuncfix->SetParName(4, "Massa");
	//fitfuncfix->FixParameter(4, Massa);

	TF1  *sigfunc = new TF1("sigfunc",BreitWigner,0 , 1.5 , 3);
	TF1  *bkgfunc = new TF1("bkgfunc",pol2,0 , 1.5 , 3);

	TCanvas *LS = new TCanvas();
	LS->SetGrid();
	BreitHisto->Draw();
	BreitHisto->Fit("fitfunc", "PQ"); //LS 
	
	
	sigfunc->FixParameter(0,fitfunc->GetParameter(3));
	sigfunc->FixParameter(1,fitfunc->GetParameter(4));
	sigfunc->FixParameter(2,fitfunc->GetParameter(5));
	sigfunc->Draw("SAME");

	
	bkgfunc->FixParameter(0,fitfunc->GetParameter(0));
	bkgfunc->FixParameter(1,fitfunc->GetParameter(1));
	bkgfunc->FixParameter(2,fitfunc->GetParameter(2));
	bkgfunc->Draw("SAME");

	TCanvas *LSF = new TCanvas();
	LSF->SetGrid();
	BreitHisto1->Draw();
	BreitHisto1->Fit("fitfuncfix", "PQ"); //LS fixed

	sigfunc->FixParameter(0,fitfunc->GetParameter(3));
	sigfunc->FixParameter(1,fitfunc->GetParameter(4));
	sigfunc->FixParameter(2,fitfunc->GetParameter(5));
	sigfunc->Draw("SAME");

	
	bkgfunc->FixParameter(0,fitfunc->GetParameter(0));
	bkgfunc->FixParameter(1,fitfunc->GetParameter(1));
	bkgfunc->FixParameter(2,NSig);
	bkgfunc->Draw("SAME");
	
	TCanvas *MLS = new TCanvas();
	MLS->SetGrid();
	BreitHisto2->Draw();
	BreitHisto2->Fit("fitfunc" , "Q"); //MLS

	sigfunc->FixParameter(0,fitfunc->GetParameter(3));
	sigfunc->FixParameter(1,fitfunc->GetParameter(4));
	sigfunc->FixParameter(2,fitfunc->GetParameter(5));
	sigfunc->Draw("SAME");

	
	bkgfunc->FixParameter(0,fitfunc->GetParameter(0));
	bkgfunc->FixParameter(1,fitfunc->GetParameter(1));
	bkgfunc->FixParameter(2,fitfunc->GetParameter(2));
	bkgfunc->Draw("SAME");

	TCanvas *MLSF = new TCanvas();
	MLSF->SetGrid();
	BreitHisto3->Draw();
	BreitHisto3->Fit("fitfuncfix" , "Q"); //MLS fixed

	sigfunc->FixParameter(0,fitfunc->GetParameter(3));
	sigfunc->FixParameter(1,fitfunc->GetParameter(4));
	sigfunc->FixParameter(2,fitfunc->GetParameter(5));
	sigfunc->Draw("SAME");

	
	bkgfunc->FixParameter(0,fitfunc->GetParameter(0));
	bkgfunc->FixParameter(1,fitfunc->GetParameter(1));
	bkgfunc->FixParameter(2,NSig);
	bkgfunc->Draw("SAME");


	TCanvas *ML = new TCanvas();
	ML->SetGrid();
	BreitHisto4->Draw();
	BreitHisto4->Fit("fitfunc", "LQ"); //ML

	sigfunc->FixParameter(0,fitfunc->GetParameter(3));
	sigfunc->FixParameter(1,fitfunc->GetParameter(4));
	sigfunc->FixParameter(2,fitfunc->GetParameter(5));
	sigfunc->Draw("SAME");

	
	bkgfunc->FixParameter(0,fitfunc->GetParameter(0));
	bkgfunc->FixParameter(1,fitfunc->GetParameter(1));
	bkgfunc->FixParameter(2,fitfunc->GetParameter(2));
	bkgfunc->Draw("SAME");
	
	
	
	

	
	
	
	
	
	gStyle->SetOptFit(1111);
	return;	
}

