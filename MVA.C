#define NSig 1e4
#define NBkg 1e6
#define bin 500
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

void MVA()
{
    Reset();	
	//Informazioni statistiche da stampare
	//gStyle->SetOptFit(1111);

    TNtuple *eve = new TNtuple("events", "events", "x:y:s"); //ebbene si, usiamo le NTuple

    TRandom3 *Rand = new TRandom3(time(0)); 

    double x, y; 
    int s;

    for (int i = 0; i<NSig; i++) //Genera il segnale
    {
        x = Rand->BreitWigner(0,1);
        y = Rand->BreitWigner(1,1);
        if(abs(x)<20 && abs(y)<20) //per ottenere un grafico leggibile dobbiamo restringere attorno a 0
            eve->Fill(x,y,1);
        else i--;
    }

    for (int i = 0; i<NBkg; i++) //Genera il background
    {
        //per generare una coppia di variabili correlate con correlazione rho 
        //ci serve prima una coppia di variabili scorrelate con varianza unitaria
        x = Rand->Gaus(0,1);
        y = Rand->Gaus(0,1);
        if(abs(x)<20 && abs(y)<20) //solito controllo
        {
            //quindi la combinazione di x
            y = y*sqrt(1-(rho*rho)) + x*rho;
            //sarà correlata con indice rho e manterrà varianza unitaria
        
            //centriamo ora la gaussiana in 4,4
            x+=4;
            y+=4;
            eve->Fill(x,y,0);
        }
        else i--;
    }

    TFile* MyFile= new TFile("MVA.root" , "recreate");
    if ( MyFile->IsOpen() ) cout<<"File opened successfully\n";
    else {
        cout<<"Error opening MVA.root\n"; 
        return;
    }
    eve->Write("MVA");

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
    TCanvas *c1 = new TCanvas();
    TCut snip = "(y<(1) || x<(0))";

    eve->Draw("y:x >> Sig2D", "s>0" && snip, "");
    eve->Draw("y:x >> Eff2D", "s>0" && !snip, "SAME");
    eve->Draw("y:x >> Bkg2D", "s<1" && !snip, "SAME");
    eve->Draw("y:x >> Pur2D", "s<1" && snip, "SAME");
    

    TH2D *Eventi = new TH2D ("Eventi" , "Eventi" , bin , -20 , 20 , bin , -20 , 20);
	Eventi->GetXaxis()->SetTitle("X");
	Eventi->GetYaxis()->SetTitle("Y");
	Eventi->SetMarkerColor(kBlue);



    TCanvas *c2 = new TCanvas();
    eve->Draw("y:x >> Eventi" , "","COLZ");


    TH1D *EventiX = new TH1D ("EventiX" , "Proiezione x" , bin , -10 , 10);
	EventiX->GetXaxis()->SetTitle("X");
	EventiX->GetYaxis()->SetTitle("Conteggi");
	//EventiX->SetMarkerColor(kBlue);

    TH1D *EventiY = new TH1D ("EventiY" , "Proiezione y" , bin , -10 , 10);
	EventiY->GetXaxis()->SetTitle("Y");
	EventiY->GetYaxis()->SetTitle("Conteggi");

    TCanvas *c3 = new TCanvas();
    c3->Divide(2 ,1);
    c3->cd(1);
    eve->Draw("x >> EventiX");
    c3->cd(2);
    eve->Draw("y >> EventiY");






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

    //MyFile->Close();

    return;
}
