//Utilizzo .x MVA.C->avvia il programma e crea i vari Canvas. Cliccare su un punto del Canvas Effi&Reiezione per avere il grafico 2D associato al semiasse (valore x)
//Se non si vuole fare lo studio sul semiasse maggiore .x MVA.C(3)->Stampa i canvas senza fare il ciclio. Con N>3 non stampa e non fa il ciclo
//.x MVA.C(1) salva i vari canvas in .png e .root NB: i canvas creati cliccando sui punti vanno salvati manualmente.
//.x MVA.C(2) salva la ntupla in un file root. Necessario se si vuole usare mycut.C per un taglio grafico.
//Taglio grafico: show tollbar-> forbice. Fai il taglio, SaveAs -> root file. Avvia mycut.C
#define NSig 1e4
#define NBkg 1e6
#define bin 500
#define rho 0.4
const static double AngBkg=47;//21.7749;//gradi
const static double AsseMin=3.05;//21.7749;//gradi
static double a=5, passo = (a-1)/25, init=1+passo;
static long superSeed=time(0);
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

string ToString(double n , int precision){
    string i=to_string((int)n);
    string d;
    double decimal=(n-(int)n)*10;
    for(int i=0; i<precision ; i++){
        d+=to_string((int)(decimal));
        decimal=(decimal-(int)decimal)*10;
    }
    return i+"."+d;
}

string Ellisse(double x0 ,double y0, double rx,double ry , double angle){
    string ang=to_string(TMath::Pi()*angle/180);
    string X1 = "(x-"+to_string(x0)+")";
    string Y1 = "(y-"+to_string(y0)+")";
    string X=X1+"*cos("+ang+")+"+Y1+"*sin("+ang+")";
    string Y=X1+"*sin("+ang+")-"+Y1+"*cos("+ang+")";
    string el="(pow("+Y+",2)/"+to_string(ry*ry)+"+pow("+X+",2)/"+to_string(rx*rx)+">1)";
    return el;
}

double ElFunction(double x0 ,double y0, double rx,double ry, double angle){
    double ang=TMath::Pi()*angle/180;
    double x1 = x0-4;
    double y1= y0-4;
    double x=x1*cos(ang)+y1*sin(ang) , y=x1*sin(ang)-y1*cos(ang);
    return y*y/(ry*ry)+x*x/(rx*rx);
}

void Click(TNtuple* ev)
{
    auto event = gPad->GetEvent();
    if(event !=11) return;
    auto select = gPad->GetSelected();
    if(!select) return;
    
    int num=0;
    double CordX=0 , CordY=0;
    double pX , pY;
    double reiezione, Sig ,purezza, efficienza, SS, SB, BS, BB;
    float_t X,Y,S;
    
    if(select->InheritsFrom("TGraph"))
    {
        TGraph *graph = (TGraph*)select;

        pX=gPad->AbsPixeltoX(gPad->GetEventX());
        pX=gPad->PadtoX(pX);
        
        pY=gPad->AbsPixeltoY(gPad->GetEventY());
        pY=gPad->PadtoY(pY);
        
        if(((pX-init)/passo)-int((pX-init)/passo)<0.5)  num=int((pX-init)/passo);
        else num=int((pX-init)/passo)+1;
        
        pX=init+passo*num;

        TH2D *Sig2Dt = new TH2D ("" , ("Ellisse Semi-Asse Maggiore: "+ToString(pX,2)).c_str() , bin , -20 , 20 , bin , -20 , 20);
        Sig2Dt->GetXaxis()->SetTitle("X");Sig2Dt->GetYaxis()->SetTitle("Y");Sig2Dt->SetMarkerColor(kGreen);Sig2Dt->SetMarkerSize(2);Sig2Dt->SetLineColor(kGreen);

        //eventi background-background -> giallo
        TH2D *Bkg2Dt = new TH2D ("" , "Bkg2D" , bin , -20 , 20 , bin , -20 , 20);
        Bkg2Dt->GetXaxis()->SetTitle("X");Bkg2Dt->GetYaxis()->SetTitle("Y");Bkg2Dt->SetMarkerColor(kYellow);Bkg2Dt->SetMarkerSize(2);Bkg2Dt->SetLineColor(kYellow);

        //eventi background-segnale -> rosso
        TH2D *Pur2Dt = new TH2D ("" , "Pur2D" , bin , -20 , 20 , bin , -20 , 20);
        Pur2Dt->GetXaxis()->SetTitle("X");Pur2Dt->GetYaxis()->SetTitle("Y");Pur2Dt->SetMarkerColor(kRed);Pur2Dt->SetMarkerSize(2);Pur2Dt->SetLineColor(kRed);

        //eventi segnale-background -> blu
        TH2D *Eff2Dt = new TH2D ("" , "Eff2D" , bin , -20 , 20 , bin , -20 , 20);
        Eff2Dt->GetXaxis()->SetTitle("X");Eff2Dt->GetYaxis()->SetTitle("Y");Eff2Dt->SetMarkerColor(kBlue);Eff2Dt->SetMarkerSize(2);Eff2Dt->SetLineColor(kBlue);

        //Proiezione x
        TH1D *AllX = new TH1D ("" , "Ellisse Proiezione X" , bin , -20 , 20);
        AllX->GetXaxis()->SetTitle("X");AllX->SetMarkerColor(kBlue);AllX->SetMarkerSize(2);AllX->SetLineColor(kBlue);

        //Proiezioney
        TH1D *AllY = new TH1D ("" , "Ellisse Proiezione Y" , bin , -20 , 20);
        AllY->GetXaxis()->SetTitle("Y");AllY->SetMarkerColor(kBlue);AllY->SetMarkerSize(2);AllY->SetLineColor(kBlue);

        ev->SetBranchAddress("x",&X);
        ev->SetBranchAddress("y",&Y);
        ev->SetBranchAddress("s",&S);

        Int_t nentries = (Int_t)ev->GetEntries();
        for (Int_t i=0;i<nentries;i++) {
            ev->GetEntry(i);
            if(S>0){
                if(ElFunction(X, Y, pX, AsseMin ,AngBkg)>1){
                    Sig2Dt->Fill(X,Y);
                    AllX->Fill(X);
                    AllY->Fill(Y);
                }
                else{
                    Eff2Dt->Fill(X,Y);
                }
            }
            else if(S<1){
                if(ElFunction(X, Y, pX, AsseMin , AngBkg)>1){
                    Pur2Dt->Fill(X,Y);
                    AllX->Fill(X);
                    AllY->Fill(Y);
                }
                else{
                    Bkg2Dt->Fill(X,Y);
                }
            }
        } 

        SS = Sig2Dt->GetEntries();
        SB = Eff2Dt->GetEntries();
        BS = Pur2Dt->GetEntries();
        BB = Bkg2Dt->GetEntries();

        purezza = 100*(SS)/(SS+BS);
        efficienza = 100*(SS)/(SS+SB);
        Sig=SS/sqrt(BS+SS);
        reiezione=100*(BB)/(BB+BS);

        TLegend* legend = new TLegend(0.1,0.7,0.4,0.9); 

        legend->SetHeader(("Ellisse p: "+ToString(purezza , 1)+"% r: "+ToString(reiezione , 1)+"%").c_str(),"C");

        legend->AddEntry(Sig2Dt , "Segnale-Segnale","lp");
        legend->AddEntry(Eff2Dt , "Segnale-Fondo","lp");
        legend->AddEntry(Bkg2Dt , "Fondo-Segnale","lp");
        legend->AddEntry(Pur2Dt , "Fondo-Fondo","lp");

        gStyle->SetOptStat(0);
        
        TCanvas* c_temp=new TCanvas();

        Sig2Dt->Draw("same");
        Eff2Dt->Draw("same");
        Bkg2Dt->Draw("same");
        Pur2Dt->Draw("same");
        legend->Draw();


        cout<<"##################################################################################\n";
        cout<<"Taglio: "<<(Ellisse(4, 4, pX, AsseMin, AngBkg)).c_str()<<endl;
        cout<<"Asse: "<<pX<<endl;

        cout<<"Segnale-Segnale: "<<SS<<endl;
        cout<<"Segnale-Fondo: "<<SB<<endl;
        cout<<"Fondo-Fondo: "<<BB<<endl;
        cout<<"Fondo-Segnale: "<<BS<<endl;

        cout <<"Purezza del Segnale: "<<purezza<<"%"<<endl;
        cout <<"Efficienza del Segnale: "<<efficienza<<"%"<<endl;
        cout <<"Reiezione del fondo: "<<reiezione<<"%"<<endl;
        cout <<"Significatività del taglio: "<<Sig<<endl;
        cout<<"##################################################################################\n";    

        TCanvas* cp = new TCanvas();
        cp->Divide(2 , 1);
        cp->cd(1);
        AllX->Draw();
        cp->cd(2);
        AllY->Draw();
        
        return;
    }
    return;
}

void MVA(int stamp=0)
{
    Reset();	
    gStyle->SetOptStat(0);
	//Informazioni statistiche da stampare
	//gStyle->SetOptFit(1111)

    TNtuple *eve = new TNtuple("events", "events", "x:y:s"); //ebbene si, usiamo le NTuple
    
    //Usare superSeed per far funzionare tutto
    TRandom3 *Rand = new TRandom3(time(0)); 

    double x, y; 
    double reiezione, Sig ,purezza, efficienza, SS, SB, BS, BB;
    int s;

    //ofstream eventsTxt;
	//eventsTxt.open("events.txt");
    for (int i = 0; i<NSig; i++) //Genera il segnale
    {
        x = Rand->BreitWigner(0,1);
        y = Rand->BreitWigner(1,1);
        if(abs(x)<20 && abs(y)<20){ //per ottenere un grafico leggibile dobbiamo restringere attorno a 0
            eve->Fill(x,y,1);
            //eventsTxt<<x<<" "<<y<<" "<<1<<endl;
        }
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
            //eventsTxt<<x<<" "<<y<<" "<<0<<endl;
        }
        else i--;
    }
    //eventsTxt.close();

    if(stamp==2){
        TFile* MyFile= new TFile("MVA.root" , "recreate");
        if ( MyFile->IsOpen() ) cout<<"File opened successfully\n";
        else {
            cout<<"Error opening MVA.root\n"; 
            return;
        }
        eve->Write("MVA");
        cout<<"File saved successfully\n";
    }

    //per poter plottare gli eventi con 4 colori diversi definisco 4 TH2D

    //###################Iperbole###########################################################################################################################
    //eventi segnale-segnale -> verde
    TH2D *Sig2D = new TH2D ("Sig2D" , "Iperbole" , bin , -20 , 20 , bin , -20 , 20);
	Sig2D->GetXaxis()->SetTitle("X");Sig2D->GetYaxis()->SetTitle("Y");Sig2D->SetMarkerColor(kGreen);Sig2D->SetMarkerSize(2);Sig2D->SetLineColor(kGreen);

    //eventi background-background -> giallo
    TH2D *Bkg2D = new TH2D ("Bkg2D" , "Bkg2D" , bin , -20 , 20 , bin , -20 , 20);
	Bkg2D->GetXaxis()->SetTitle("X");Bkg2D->GetYaxis()->SetTitle("Y");Bkg2D->SetMarkerColor(kYellow);Bkg2D->SetMarkerSize(2);Bkg2D->SetLineColor(kYellow);

    //eventi background-segnale -> rosso
    TH2D *Pur2D = new TH2D ("Pur2D" , "Pur2D" , bin , -20 , 20 , bin , -20 , 20);
	Pur2D->GetXaxis()->SetTitle("X");Pur2D->GetYaxis()->SetTitle("Y");Pur2D->SetMarkerColor(kRed);Pur2D->SetMarkerSize(2);Pur2D->SetLineColor(kRed);

    //eventi segnale-background -> blu
    TH2D *Eff2D = new TH2D ("Eff2D" , "Eff2D" , bin , -20 , 20 , bin , -20 , 20);
	Eff2D->GetXaxis()->SetTitle("X");Eff2D->GetYaxis()->SetTitle("Y");Eff2D->SetMarkerColor(kBlue);Eff2D->SetMarkerSize(2);Eff2D->SetLineColor(kBlue);
    
    //########Snip######################################################################################################################################
    //eventi segnale-segnale -> verde
    TH2D *Sig2DS = new TH2D ("Sig2DS" , "Taglio Semplice" , bin , -20 , 20 , bin , -20 , 20);
	Sig2DS->GetXaxis()->SetTitle("X");Sig2DS->GetYaxis()->SetTitle("Y");Sig2DS->SetMarkerColor(kGreen);Sig2DS->SetMarkerSize(2);Sig2DS->SetLineColor(kGreen);

    //eventi background-background -> giallo
    TH2D *Bkg2DS = new TH2D ("Bkg2DS" , "Bkg2DS" , bin , -20 , 20 , bin , -20 , 20);
	Bkg2DS->GetXaxis()->SetTitle("X");Bkg2DS->GetYaxis()->SetTitle("Y");Bkg2DS->SetMarkerColor(kYellow);Bkg2DS->SetMarkerSize(2);Bkg2DS->SetLineColor(kYellow);

    //eventi background-segnale -> rosso
    TH2D *Pur2DS = new TH2D ("Pur2DS" , "Pur2DS" , bin , -20 , 20 , bin , -20 , 20);
	Pur2DS->GetXaxis()->SetTitle("X");Pur2DS->GetYaxis()->SetTitle("Y");Pur2DS->SetMarkerColor(kRed);Pur2DS->SetMarkerSize(2);Pur2DS->SetLineColor(kRed);

    //eventi segnale-background -> blu
    TH2D *Eff2DS = new TH2D ("Eff2DS" , "Eff2DS" , bin , -20 , 20 , bin , -20 , 20);
	Eff2DS->GetXaxis()->SetTitle("X");Eff2DS->GetYaxis()->SetTitle("Y");Eff2DS->SetMarkerColor(kBlue);Eff2DS->SetMarkerSize(2);Eff2DS->SetLineColor(kBlue);

    //#################################################################################################################################################
    //variabili per i cut
    string basiccut="(y<0 || x<1)";
    string ipercut="(y<1/(x)+1 || y<0 || x<0)";

    //Cut testati
    TCut ip = ipercut.c_str();
    TCut snip =basiccut.c_str();//(Ellisse(4, 4 , 2, 3)).c_str();

    //disegno sullo stesso canvas gli eventi, divisi in modo opportuno rispetto al taglio
    TCanvas *c1 = new TCanvas();

    //Primo taglio
    eve->Draw("y:x >> Sig2D", "s>0" && ip , "");
    eve->Draw("y:x >> Eff2D", "s>0" && !ip , "same");
    eve->Draw("y:x >> Bkg2D", "s<1" && !ip , "same");
    eve->Draw("y:x >> Pur2D", "s<1" && ip , "same");

    SS = Sig2D->GetEntries();
    SB = Eff2D->GetEntries();
    BS = Pur2D->GetEntries();
    BB = Bkg2D->GetEntries();

    purezza = 100*(SS)/(SS+BS);
    efficienza = 100*(SS)/(SS+SB);
    Sig=SS/sqrt(BS+SS);
    reiezione=100*(BB)/(BB+BS);

    TLegend* legendIP = new TLegend(0.1,0.7,0.4,0.9); 
    
    legendIP->SetHeader(("Iperbole p: "+ToString(purezza , 1)+"% r: "+ToString(reiezione , 1)+"%").c_str(),"C");

    legendIP->AddEntry(Sig2D , "Segnale-Segnale","lp");legendIP->AddEntry(Eff2D , "Segnale-Fondo","lp");
    legendIP->AddEntry(Bkg2D , "Fondo-Segnale","lp");legendIP->AddEntry(Pur2D , "Fondo-Fondo","lp");
    legendIP->Draw();

    if (stamp==1 ||stamp==3){
        int Cartella1= system("mkdir -p MVA");
        c1->SaveAs("MVA/iperbole.png");
        c1->SaveAs("MVA/iperbole.root");
    }

    cout<<"##################################################################################\n";
    cout<<"Taglio: "<<ipercut.c_str()<<endl;

    cout<<"Segnale-Segnale: "<<SS<<endl;
    cout<<"Segnale-Fondo: "<<SB<<endl;
    cout<<"Fondo-Fondo: "<<BB<<endl;
    cout<<"Fondo-Segnale: "<<BS<<endl;

    cout <<"Purezza del Segnale Iperbole: "<<purezza<<"%"<<endl;
    cout <<"Efficienza del Segnale Iperbole: "<<efficienza<<"%"<<endl;
    cout <<"Reiezione del fondo Iperbole: "<<reiezione<<"%"<<endl;
    cout <<"Significatività del taglio Iperbole: "<<Sig<<endl;
    cout<<"##################################################################################\n\n";

    //Secondo taglio
    TCanvas* cs = new TCanvas();

    eve->Draw("y:x >> Sig2DS", "s>0" && snip, "");
    eve->Draw("y:x >> Eff2DS", "s>0" && !snip, "SAME");
    eve->Draw("y:x >> Bkg2DS", "s<1" && !snip, "SAME");
    eve->Draw("y:x >> Pur2DS", "s<1" && snip, "SAME");

    //calcolo purezza ed efficienza del taglio

    SS = Sig2DS->GetEntries();
    SB = Eff2DS->GetEntries();
    BS = Pur2DS->GetEntries();
    BB = Bkg2DS->GetEntries();

    purezza = 100*(SS)/(SS+BS);
    efficienza = 100*(SS)/(SS+SB);
    Sig=SS/sqrt(BS+SS);
    reiezione=100*(BB)/(BB+BS);

    TLegend* legend = new TLegend(0.1,0.7,0.4,0.9); 

    legend->SetHeader(("Taglio Semplice p: "+ToString(purezza , 1)+"% r: "+ToString(reiezione , 1)+"%").c_str(),"C");

    legend->AddEntry(Sig2DS , "Segnale-Segnale","lp");legend->AddEntry(Eff2DS , "Segnale-Fondo","lp");
    legend->AddEntry(Bkg2DS , "Fondo-Segnale","lp");legend->AddEntry(Pur2DS , "Fondo-Fondo","lp");
    legend->Draw();

    if (stamp==1 ||stamp==3){
        //int Cartella= system("mkdir -p MVA");
        c1->SaveAs("MVA/snip.png");
        c1->SaveAs("MVA/snip.root");
    }

    cout<<"##################################################################################\n";
    cout<<"Taglio: "<<basiccut.c_str()<<endl;

    cout<<"Segnale-Segnale: "<<SS<<endl;
    cout<<"Segnale-Fondo: "<<SB<<endl;
    cout<<"Fondo-Fondo: "<<BB<<endl;
    cout<<"Fondo-Segnale: "<<BS<<endl;

    cout <<"Purezza del Segnale Snip: "<<purezza<<"%"<<endl;
    cout <<"Efficienza del Segnale Snip: "<<efficienza<<"%"<<endl;
    cout <<"Reiezione del fondo Snip: "<<reiezione<<"%"<<endl;
    cout <<"Significatività del taglio Snip: "<<Sig<<endl;
    cout<<"##################################################################################\n";



//###################################### Ellisse varibile come taglio: Test ######################################################
    if(stamp==0 || stamp==1){
        TGraph* segnale = new TGraph();
        segnale->SetLineColor(kRed);
        segnale->SetMarkerStyle(kStar);

        TGraph* fondo = new TGraph();
        fondo->SetMarkerStyle(kPlus);
        fondo->SetLineColor(kBlue);
        TCut ellisse;

        TMultiGraph* mg = new TMultiGraph();

        TH2D *Sig2DT = new TH2D ("Sig2DT" , "Plot Run" , bin , -20 , 20 , bin , -20 , 20);
        Sig2DT->GetXaxis()->SetTitle("X");Sig2DT->GetYaxis()->SetTitle("Y");Sig2DT->SetMarkerColor(kGreen);

        //eventi background-background -> giallo
        TH2D *Bkg2DT = new TH2D ("Bkg2DT" , "Bkg2D" , bin , -20 , 20 , bin , -20 , 20);
        Bkg2DT->GetXaxis()->SetTitle("X");Bkg2DT->GetYaxis()->SetTitle("Y");Bkg2DT->SetMarkerColor(kYellow);

        //eventi background-segnale -> rosso
        TH2D *Pur2DT = new TH2D ("Pur2DT" , "Pur2D" , bin , -20 , 20 , bin , -20 , 20);
        Pur2DT->GetXaxis()->SetTitle("X");Pur2DT->GetYaxis()->SetTitle("Y");Pur2DT->SetMarkerColor(kRed);

        //eventi segnale-background -> blu
        TH2D *Eff2DT = new TH2D ("Eff2DT" , "Eff2D" , bin , -20 , 20 , bin , -20 , 20);
        Eff2DT->GetXaxis()->SetTitle("X");Eff2DT->GetYaxis()->SetTitle("Y");Eff2DT->SetMarkerColor(kBlue);

        int j=0;
        for(double i =1; i<=a; i=i+passo){

            ellisse = (Ellisse(4, 4, i, AsseMin , AngBkg)).c_str();

            //cout<< (Ellisse(4, 4, a, a*(0.88))).c_str()<<endl;

            eve->Draw("y:x >> Sig2DT", "s>0" && ellisse, "goff");
            eve->Draw("y:x >> Eff2DT", "s>0" && !ellisse, "goff");
            eve->Draw("y:x >> Bkg2DT", "s<1" && !ellisse, "goff");
            eve->Draw("y:x >> Pur2DT", "s<1" && ellisse, "goff");

            SS = Sig2DT->GetEntries();
            SB = Eff2DT->GetEntries();
            BS = Pur2DT->GetEntries();
            BB = Bkg2DT->GetEntries();

            purezza = 100*(SS)/(SS+BS);
            efficienza = 100*(SS)/(SS+SB);
            Sig=SS/sqrt(BS+SS);
            reiezione=100*(BB)/(BB+BS);

            cout<<j+1<<") "<<"Efficienza= "<<efficienza<<" , Reiezione="<<reiezione<<" , Asse: "<<i<<" , Significatività: "<<Sig<<endl;

            segnale->SetPoint(j,i , efficienza/100);
            fondo->SetPoint(j,i , reiezione/100);
            j++;

        }
        TLegend* legendg = new TLegend(); 
        legendg->AddEntry(segnale , "Efficienza Segnale","lp");legendg->AddEntry(fondo , "Reiezione Fondo","lp");

        TCanvas* c_prova = new TCanvas();
        mg->Add(segnale);
        mg->Add(fondo);
        mg->SetTitle("Roc; SemiAsse; Efficienza");

        c_prova->SetGrid();
        mg->Draw("APC");
        legendg->Draw();
        gPad->AddExec("ex" , "Click(events)");
        if(stamp==1){
            c_prova->SaveAs("MVA/Roc.png");
            // c_prova->SaveAs("MVA/Roc.root");//La funzione di click fa macello se salvi il .root
        }
    }

//################################################################################################################################
    TH2D *Eventi = new TH2D ("Eventi" , "Eventi" , bin , -20 , 20 , bin , -20 , 20);
	Eventi->GetXaxis()->SetTitle("X");
	Eventi->GetYaxis()->SetTitle("Y");
	Eventi->SetMarkerColor(kBlue);

    TH1D *EventiX = new TH1D ("EventiX" , "Proiezione x" , bin , -10 , 10);
	EventiX->GetXaxis()->SetTitle("X");
	EventiX->GetYaxis()->SetTitle("Conteggi");
	//EventiX->SetMarkerColor(kBlue);

    TH1D *EventiY = new TH1D ("EventiY" , "Proiezione y" , bin , -10 , 10);
	EventiY->GetXaxis()->SetTitle("Y");
	EventiY->GetYaxis()->SetTitle("Conteggi");

    TCanvas *c2 = new TCanvas();
    eve->Draw("y:x >> Eventi" , "","COLZ");

    TCanvas *c3 = new TCanvas();
    c3->Divide(2 ,1);
    c3->cd(1);
    eve->Draw("x >> EventiX");
    c3->cd(2);
    eve->Draw("y >> EventiY");

    if (stamp==1 ||stamp==3){
        c2->SaveAs("MVA/AllEv.png");
        c2->SaveAs("MVA/AllEv.root");
        c3->SaveAs("MVA/Proizione.png");
        c3->SaveAs("MVA/Proiezione.root");

    }
    //MyFile->Close();

    return;
}
