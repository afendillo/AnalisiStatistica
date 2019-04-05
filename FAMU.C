#define Pi TMath::Pi()
#define NFibre 32
#define spessore 0.3 //cm
#define lunghezza 10 //cm
#define FAC 0.04
#define Q_e 1.60217662e-19
#define L 210 //cm

using namespace std;
struct posizione3D{double x , y , theta, phi;int piano, fibra, flag=0;};
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
    angle[0]=(rdm->Rndm()+1)*Pi/2;//theta
    double hit=rdm->Rndm();
    while(hit>pow(cos(angle[0]), 2)){
        hit=rdm->Rndm();
        angle[0]=(rdm->Rndm()+1)*Pi/2;
    }
    angle[1]=rdm->Rndm()*2*Pi;//phi
    
    return angle;
}

posizione3D CoordPoint1P(TRandom3* rdm)
{
    //origine estremo del piano quadrato-->prima fibra = pos 0
    posizione3D pos;
    pos.x=rdm->Rndm()*(NFibre)*spessore;
    pos.y=rdm->Rndm()*lunghezza;
    double *posTP=Angle(rdm);
    pos.theta=posTP[0];
    pos.phi=posTP[1];
    if (pos.theta!=0)
    {
        if((pos.x/spessore)-(int)(pos.x/spessore)==0 && pos.phi>Pi )
        {
            pos.fibra=ceil(pos.x/spessore)-1;
        }
        else pos.fibra=ceil(pos.x/spessore);//fibra colpita
        pos.piano=1;
    }
    
    return pos;
}

posizione3D Proiezione(posizione3D pos)
{
    posizione3D NewPos;
    NewPos=pos;
    if (pos.theta!=0 && pos.theta!=Pi/2)
    {
        NewPos.x=pos.x-1*sin(NewPos.phi)*spessore/tan(NewPos.theta);//Asse Z positivo verso l'alto, i muoni vanno verso il basso
        NewPos.y=pos.y-1*cos(NewPos.phi)*spessore/tan(NewPos.theta);
        if((NewPos.x/spessore)-(int)(NewPos.x/spessore)==0 && NewPos.phi>Pi)
        {
            NewPos.fibra=ceil(pos.x/spessore)-1;
        }
        else NewPos.fibra=ceil(NewPos.x/spessore);//fibra colpita
        
        NewPos.piano=1;
    }
    else if(pos.theta==Pi/2) NewPos=pos;

    if(NewPos.x>9.6 || NewPos.x<0 || NewPos.y<0 || NewPos.y>10) NewPos.flag=1;
    return NewPos;
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
    double Conv = 180/Pi;

    TRandom3* rand= new TRandom3(time(0));
    posizione3D posPiano1Up, posPiano1D , posPiano2D, posPiano2Up;
    for(int i=0; i<10;i++){
        posPiano1Up = CoordPoint1P(rand);
        posPiano1D= Proiezione(posPiano1Up);
        posPiano2Up =posPiano1D;
        posPiano2Up.x=posPiano1D.y;
        posPiano2Up.y=posPiano1D.x;
        //double *p=CoordinateFibra12Theta(posPiano1Up);
        cout<<"################################################################################################################## \n";
        cout<<"\n PIANO "<<posPiano1Up.piano<<" Up \n";
        cout<<"\n X : "<<posPiano1Up.x<<" , Y : "<<posPiano1Up.y<<" Angolo theta: "<< Conv*posPiano1Up.theta<<"° , Angolo phi: "<<Conv*posPiano1Up.phi<<"°"<<endl;
        cout<<"Fibra numero: "<<posPiano1Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano1Up.y<<"\n";

        cout<<"\n PIANO "<<posPiano1D.piano<<" D \n";
        cout<<"\n X : "<<posPiano1D.x<<" , Y : "<<posPiano1D.y<<" Angolo theta: "<<Conv*posPiano1D.theta<<"° , Angolo phi: "<<Conv*posPiano1D.phi<<"° \n";
        if (posPiano1D.flag==0) cout<<"Fibra numero: "<<posPiano1D.fibra<<" , distanza da PM: "<<lunghezza-posPiano1D.y<<"\n";
        else cout<<"Fuori accettanza \n";
        cout<<"------------------------------------------------ \n";
        // posPiano2D=Proiezione(posPiano2Up);
        // cout<<"\n PIANO "<<posPiano2Up.piano<<" D \n";
        // cout<<"\n X : "<<posPiano2Up.x<<" , Y : "<<posPiano2Up.y<<" Angolo theta: "<< posPiano2Up.theta<<" , Angolo phi: "<<posPiano2Up.phi<<"\n";
        // if (posPiano2Up.flag==0) cout<<"Fibra numero: "<<posPiano2Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano2Up.y<<"\n";
        // else cout<<"Fuori accettanza \n";


        cout<<"################################################################################################################## \n";
        //delete [] p;
    }
     return;
}