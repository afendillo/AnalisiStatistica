#define NSig 10000
#define NBkg 1000000
#define bin 1000
#define rho 0.4
using namespace std;
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

void mycut(char stamp='Q')
{
	Reset();
	//Informazioni statistiche da stampare
	//gStyle->SetOptFit(1111);

    double PlotRun;

    TFile* MyFile = new TFile("MVA.root");

    if ( MyFile->IsOpen() ) cout<<"File opened successfully\n";
    else return;

    auto eve = (TNtuple *)((MyFile->GetKey("MVA"))->ReadObj());

    TFile* MyCut = new TFile("mycut.root");
    
    if ( MyCut->IsOpen() ) cout<<"Cut opened successfully\n";
    else return;

    auto cutg = (TCutG *)((MyCut->GetKey("CUTG"))->ReadObj());
    cutg->SetVarX("x");
    cutg->SetVarY("y");

    //per poter plottare gli eventi con 4 colori diversi definisco 4 TH2D
    
    //eventi segnale-segnale -> verde
    TH2D *Sig2D = new TH2D ("Sig2D" , "Plot Run" , bin , -20 , 20 , bin , -20 , 20);
	Sig2D->GetXaxis()->SetTitle("X");
	Sig2D->GetYaxis()->SetTitle("Y");
	Sig2D->SetMarkerColor(kGreen);

    //eventi background-background -> giallo
    TH2D *Bkg2D = new TH2D ("Bkg2D" , "Bkg2D" , bin , -20 , 20 , bin , -20 , 20);
	Bkg2D->GetXaxis()->SetTitle("X");
	Bkg2D->GetYaxis()->SetTitle("Y");
	Bkg2D->SetMarkerColor(kYellow);

    //eventi background-segnale -> rosso
    TH2D *Pur2D = new TH2D ("Pur2D" , "Pur2D" , bin , -20 , 20 , bin , -20 , 20);
	Pur2D->GetXaxis()->SetTitle("X");
	Pur2D->GetYaxis()->SetTitle("Y");
	Pur2D->SetMarkerColor(kRed);

    //eventi segnale-background -> blu
    TH2D *Eff2D = new TH2D ("Eff2D" , "Eff2D" , bin , -20 , 20 , bin , -20 , 20);
	Eff2D->GetXaxis()->SetTitle("X");
	Eff2D->GetYaxis()->SetTitle("Y");
	Eff2D->SetMarkerColor(kBlue);

    //disegno sullo stesso canvas gli eventi, divisi in modo opportuno rispetto al taglio
    TCanvas *c2 = new TCanvas();

    eve->Draw("y:x >> Sig2D","s>0 && !CUTG", "");
    eve->Draw("y:x >> Eff2D", "s>0 && CUTG", "SAME");
    eve->Draw("y:x >> Bkg2D", "s<1 && CUTG", "SAME");
    eve->Draw("y:x >> Pur2D", "s<1 && !CUTG", "SAME");
    
    //calcolo purezza ed efficienza del taglio

    double purezza, efficienza, SS, SB, BS, BB;

    SS = Sig2D->GetEntries();
    SB = Eff2D->GetEntries();
    BS = Pur2D->GetEntries();
    BB = Bkg2D->GetEntries();

    purezza = 100*(SS)/(SS+BS);
    efficienza = 100*(SS)/(SS+SB);

    cout <<"Purezza del taglio: "<<purezza<<"%"<<endl;
    cout <<"Efficienza del taglio: "<<efficienza<<"%"<<endl;
 
    return;
}