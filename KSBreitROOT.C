#define Pi TMath::Pi() 
#define Massa 91.187
#define Gamma 2.495
#define Range 4
#define NRand 1e6
#define bin 100



using namespace std;

//bad rando
static const double C = 65539;
static const double RANDU_MAX = pow(2 , 31);
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
double Randu()
{
    static double x = 1;
    x = fmod(C*x, RANDU_MAX);
    return x/RANDU_MAX;
}
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
	
	double y = 1/p[0]*(2*Gamma*Range)*NRand/bin*(Gamma/2)/(pow(x[0]-Massa,2)+Gamma*Gamma/4);
	return y;
}

//Constant PDFunction a-b interval
double cFunc(double x) {
		return ROOT::Math::uniform_pdf  (x,0,1,0);
							
}	

void KSBreitROOT(char stamp='Q')
{	
	Reset();
	
	//Informazioni statistiche da stampare
	gStyle->SetOptFit(1111);
//----------------------------------------------------------------------------------------------------------------------	
//------------------------------Inverse function method for Normal Distribution-----------------------------------------	
//----------------------------------------------------------------------------------------------------------------------

//Paragone fra Gauss_InverseMethod e Random3->Gaus	
	// int Media=10, Sigma=1;
	
	// TRandom3 *Gauss = new TRandom3((unsigned)time( NULL ));
	// TRandom3 *Gauss2 = new TRandom3((unsigned)time( NULL ));

	// TH1F *gauss = new TH1F ("gauss" , "gauss" , 100 , 0, Media*2);
	// gauss->SetFillColorAlpha(kBlue, 0.30);
	
	// TH1F *gauss2 = new TH1F ("gauss2" , "gauss2" , 100 , 0, Media*2);
	// gauss2->SetFillColorAlpha(kRed, 0.30);
	
	// TH1F *gaussC = new TH1F ("gaussC" , "gaussC" , 100 , 0, Media*2);
	// gaussC->SetFillColorAlpha(kRed, 0.30);
	
	// TF1 *f1 = new TF1("f1", "[0]*[2]*exp(-0.5*(pow((x-[1])/[2] , 2)))", 0, Media*2);
	// f1->SetParameter(0 , sqrt(TMath::TwoPi()));
	// f1->SetParameter(1 , Media);
	// f1->SetParameter(2 , Sigma);
	
	// TF1 *f2 = new TF1("f2", "[0]+[1]*sqrt(2)*myErfInv(2*x-1)", Media, Media*2);
	// f2->SetParameter(0 , Media);
	// f2->SetParameter(1 , Sigma);
	
	// for (int i = 0 ; i<1000; i++)
	// {
		// gauss->Fill(Gauss->Gaus(Media,Sigma));
		// gauss2->Fill(f2->Eval(Gauss2->Rndm()));
		// gaussC->Fill(f2->Eval((float)rand()/RAND_MAX));
		
	// }

	// TCanvas *c1 = new TCanvas();
	// c1->SetGrid();
	// gauss->Draw();
	// gauss->Fit("f1" , "Q");
	
	// TCanvas *c2 = new TCanvas();
	// c2->SetGrid();
	// gauss2->Draw();
	// gauss2->Fit("f1" , "Q");
	
	// TCanvas *c3 = new TCanvas();
	// c3->SetGrid();
	// gaussC->Draw();
	// gaussC->Fit("f1" , "Q");
	
//----------------------------------------------------------------------------------------------------------------------	
//--------------------------------- Randomness Test for TRandom3 -------------------------------------------------------	
//----------------------------------------------------------------------------------------------------------------------		
	double Max=Pi/2 , Min=-Pi/2;
	double random=0. ,randomX=0. , randomY=0. , randomZ=0.;	
	double Media=0., Media2=0.,Varianza=0.;
	
	//Se il seed non è settato verranno generati gli stessi RN ad ogni avvio
	TRandom3 *RandGenerator = new TRandom3(time(0));
	
	TH1D *Test1D = new TH1D ("Test1D" , "Test1D" , bin , 0 , 1);
	Test1D->GetXaxis()->SetTitle("X");
	Test1D->GetYaxis()->SetTitle("Conteggi");
	Test1D->SetFillColorAlpha(kRed, 0.30);
	
	TH2D *Test2D = new TH2D ("Test2D" , "Test2D" , bin , 0 , 1 , bin , 0 , 1);
	Test2D->GetXaxis()->SetTitle("X");
	Test2D->GetYaxis()->SetTitle("Conteggi");
	Test2D->SetMarkerColor(kRed);
	
	TH3D *Test3D = new TH3D ("Test3D" , "Test3D" , bin , 0 , 1 , bin , 0 , 1 , bin , -1 , 2);
	Test3D->GetXaxis()->SetTitle("X");
	Test3D->GetYaxis()->SetTitle("Conteggi");
	Test3D->SetMarkerColor(kRed);
	
    TGraph* graph2D = new TGraph();
    graph2D->SetTitle("Graph2D");
    graph2D->GetXaxis()->SetTitle("X");
	graph2D->GetYaxis()->SetTitle("Conteggi");
	graph2D->SetMarkerColor(kRed);
    
	for (int i = 0 ; i<NRand; i++)
	{
		randomX=RandGenerator->Rndm();
		randomY=RandGenerator->Rndm();
		randomZ=RandGenerator->Rndm();
		Test1D->Fill(randomX);	
		Test2D->Fill(randomX , randomY);	
		Test3D->Fill(randomX , randomY , randomZ);
        graph2D->SetPoint(i , randomX , randomY );
		
		Media+=randomX;
		Media2+=randomX*randomX;		
		
	}
	
	Media=Media/NRand;
	Varianza=sqrt((Media2/NRand-Media*Media));
	
	//cout<<(2<3)*(NRand==NRand)<<endl;
	
	TCanvas *CasTest1D = new TCanvas();
	CasTest1D->SetGrid();
	Test1D->Draw();
	
	TCanvas *CasTest2D = new TCanvas();
	CasTest2D->SetGrid();
	Test2D->Draw("E");
	
	TCanvas *CasTest3D = new TCanvas();
	CasTest3D->SetGrid();
	Test3D->Draw();
    
    TCanvas *TestGraph2D = new TCanvas();
	TestGraph2D->SetGrid();
	graph2D->Draw();
//----------------------------------------------------------------------------------------------------------------------	
//--------------------------------- Kolmogorov-Smirnov Test for Uniform Distribution -----------------------------------	
//----------------------------------------------------------------------------------------------------------------------	
	double *sample = new double[NRand];
	for (int i = 0 ; i<NRand; i++)
	{	
		//sample[i] = Randu();
		sample [i] = RandGenerator->Rndm();
	}
	ROOT::Math::Functor1D f(&cFunc);
	ROOT::Math::GoFTest *goftest = new ROOT::Math::GoFTest(NRand, sample, f, ROOT::Math::GoFTest::kPDF, 0,1);

	cout<<"Kolmogorov-Smirnov Test: "<<goftest->KolmogorovSmirnovTest()<<endl;
//----------------------------------------------------------------------------------------------------------------------	
//--------------------------------- Inverse Function Method for Breit-Wigner -------------------------------------------	
//----------------------------------------------------------------------------------------------------------------------	
	
	TF1 *breit = new TF1("breit",BreitWigner, Massa-Gamma*Range, Massa+Gamma*Range , 1);
	breit->SetParameter(0, 1);
	
	TH1D *BreitHisto = new TH1D ("BreitHisto" , "BreitHisto" , bin ,Massa-Gamma*Range, Massa+Gamma*Range);
	BreitHisto->GetXaxis()->SetTitle("Massa [GeV]");
	BreitHisto->GetYaxis()->SetTitle("Conteggi");
	BreitHisto->SetFillColorAlpha(kGreen, 0.30);
	
	for (int i = 0 ; i<NRand; i++)
	{
		random=Massa+Gamma/2*tan(RandGenerator->Rndm()*(Max-Min)+Min-TMath::ATan(2*Massa/Gamma));
		BreitHisto->Fill(random);	
	}
	
	TCanvas *test = new TCanvas();
	test->SetGrid();

	BreitHisto->Draw();

	BreitHisto->Fit("breit" , "Q");
	gStyle->SetOptFit(1111);
	return;	
}

