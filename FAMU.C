#define Pi TMath::Pi()
#define NFibre 32
#define spessore 0.3 //cm
#define lunghezza 10 //cm
#define FAC 0.04
#define Q_e 1.60217662e-19
#define L 210 //cm

using namespace std;
struct posizione3D{double x , y , theta, phi;int piano, fibra;};
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

double* Angle(TRandom3* rdm)
{
    double * angle= new double[2];
    angle[0]=rdm->Rndm()*Pi;//theta
    double hit=rdm->Rndm();
    while(hit>pow(cos(angle[0]), 2)){
        hit=rdm->Rndm();
        angle[0]=rdm->Rndm()*Pi;
    }
    angle[1]=rdm->Rndm()*2*Pi;//phi
    
    return angle;
}

posizione3D CoordPoint(TRandom3* rdm)
{
    //origine estremo del piano quadrato-->prima fibra = pos 0
    posizione3D pos;
    pos.x=rdm->Rndm()*(NFibre)*spessore;
    pos.y=rdm->Rndm()*lunghezza;
    double *posTP=Angle(rdm);
    pos.theta=posTP[0];
    pos.phi=posTP[1];
    if (pos.theta!=0 && pos.theta!=Pi)
    {
        pos.fibra=round(abs(pos.x-0.15)/0.3)+1;//fibra colpita
    }
    pos.piano=1;
    return pos;
}


// double* CoordinateFibra12Theta(posizione3D pos)
// {
//     int NF;
//     double *Coord= new double[4];
//     //Hit per il primo piano
//     if (pos.theta!=0 && pos.theta!=Pi && pos.theta!=Pi/2)
//     {
//         Coord[0]=round(abs(pos.x-0.15)/0.3)+1;//fibra colpita
//     }
//     Coord[1]=lunghezza-pos.y;//posizione y della fibra colpita
//     Coord[2]=pos.theta;
//     Coord[4]=pos.phi;
//     return Coord;
// }
void FAMU(){

    Reset();
    double I0=1;//cm-2s-1
    double G=1e6;//Gain
    double QE=0.5;//efficienza di conversione

    TRandom3* rand= new TRandom3(time(0));
    posizione3D posPiano1Up;
    for(int i=0; i<1 ; i++){
        posPiano1Up = CoordPoint(rand);
        //double *p=CoordinateFibra12Theta(posPiano1Up);
        cout<<"\nX : "<<posPiano1Up.x<<" , Y : "<<posPiano1Up.y<<" Angolo theta: "<< posPiano1Up.theta<<" , Angolo phi: "<<posPiano1Up.phi<<endl;
        cout<<"Fibra piano "<<posPiano1Up.piano<<" numero: "<<posPiano1Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano1Up.y<<endl;
        //delete [] p;
    }


     return;
}