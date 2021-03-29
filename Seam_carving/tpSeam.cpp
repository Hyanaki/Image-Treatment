#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "imagesRW.h"
#include <string.h>
#include <iostream>

// Fonction permettant de g�rer les effets de miroir sur les bords de l'image
void miroir(int &x, int &y, const unsigned char *I, int sizeX, int sizeY) {

	// Si la valeur de x est inf�rieur � 0 et est donc en dehors (au dessus visuellement) de l'image on prend 
	// la valeur r�cup�r�e par effet de miroir par rapport au bord du dessus de l'image
	if (x<0) {
		x = -x - 1; // on met la valeur n�gatif -> positive puis on soustrait 1 car la valeur 0 fait partie de l'image (bord)					
	}
	// Si la valeur de y est inf�rieur � 0 et est donc en dehors (� gauche visuellement) de l'image on prend 
	// la valeur r�cup�r�e par effet de miroir par rapport au bord gauche de l'image
	if (y<0) {
		y = -y - 1; // on met la valeur n�gatif -> positive puis on soustrait 1 car la valeur 0 fait partie de l'image (bord)
	}
	// Si la valeur de x est sup�rieur ou �gale � la hauteur de l'image et est donc en dehors (en dessous visuellement)
	// de l'image on prend la valeur r�cup�r�e par effet de miroir par rapport au bord du dessous de l'image
	if (x >= sizeY) {
		x = sizeY - 1 - (x - sizeY); // La valeur au bord moins la diff�rence entre le bord et la valeur de x
	}
	// Si la valeur de y est sup�rieur ou �gale � la largeur de l'image et est donc en dehors (� droite visuellement)
	// de l'image on prend la valeur r�cup�r�e par effet de miroir par rapport au bord droite de l'image
	if (y >= sizeX) {
		y = sizeX - 1 - (y - sizeX); // La valeur au bord moins la diff�rence entre le bord et la valeur de y
	}
}

// Calcul du gradient de imgIn  (Attention aux bords)
void IMgradient_mask (unsigned char *imgIn, unsigned char *imgOut, int sizeX, int sizeY, unsigned char *mask)
{
	// Initialisation de 2 tableau qui contiendront le gradient de l'image selon l'axe x et y
	double gradx[sizeX*sizeY];
	double grady[sizeY*sizeX];

	// Parcours de tous les pixels de l'image
	for (int i = 0; i<sizeY; i++) {
		for (int j = 0; j<sizeX; j++) {
			// Si le pixel voisin selon y est en dehors de l'image alors on assigne la valeur du pixel dans le gradient en y � la valeur du pixel de l'image
			// autrement on fait la diff�rence entre le pixel et son voisin afin d'avoir des valeurs �lev�es aux contours dans l'image
			if (i + 1>=sizeY) gradx[j + sizeX*i] = imgIn[j + sizeX*i];
			else { gradx[j + sizeX*i] = imgIn[j + sizeX*(i + 1)] - imgIn[j + sizeX*i]; }

			// Si le pixel voisin selon x est en dehors de l'image alors on assigne la valeur du pixel dans le gradient en x � la valeur du pixel de l'image
			// autrement on fait la diff�rence entre le pixel et son voisin afin d'avoir des valeurs �lev�es aux contours dans l'image
			if (j + 1>=sizeX) grady[j + sizeX*i] = imgIn[j + sizeX*i];
			else { grady[j + sizeX*i] = imgIn[j + 1 + sizeX*i] - imgIn[j + sizeX*i]; }
		}
	}

	// Afin d'obtenir le gradient de l'image, il reste plus qu'a faire la norme de l'addition des 2 matrices de gradient selon les diff�rents axes
	for (int i = 0; i<sizeY; i++) {
		for (int j = 0; j<sizeX; j++) {
			imgOut[j + sizeX*i] = sqrt(pow(gradx[j + sizeX*i], 2) + pow(grady[j + sizeX*i], 2));
		}
	}

	//MASQUAGE
	for(int i=0; i<sizeY;i++){
		for(int j=0; j<sizeX;j++){
			// Lorsque la valeur de masque fais partie de la forme on met l'�nergie du gradient � 0 afin de facilit� le passage de la coupure � cette endroit
			if(mask[i*sizeX + j]==0)
			{
				imgOut[i*sizeX + j] =0;
				// imgOut[i*sizeX + j] =1000; // Cas d'un masquage permettant de garder l'element masqu�
			}
		}
	}
}

// Calcul du gradient de imgIn  (Attention aux bords)
void IMgradient (unsigned char *imgIn, unsigned char *imgOut, int sizeX, int sizeY)
{
	// Initialisation de 2 tableau qui contiendront le gradient de l'image selon l'axe x et y
	double gradx[sizeX*sizeY];
	double grady[sizeY*sizeX];
	
	for(int i=0; i<sizeY; i++){
		for (int j=0; j<sizeX; j++){
			// Si le pixel voisin selon y est en dehors de l'image alors on assigne la valeur du pixel dans le gradient en y � la valeur du pixel de l'image
			// autrement on fait la diff�rence entre le pixel et son voisin afin d'avoir des valeurs �lev�es aux contours dans l'image
			if (i + 1 >= sizeY) gradx[j + sizeX*i] = imgIn[j + sizeX*i];
			else { gradx[j + sizeX*i] = imgIn[j + sizeX*(i + 1)] - imgIn[j + sizeX*i]; }

			// Si le pixel voisin selon x est en dehors de l'image alors on assigne la valeur du pixel dans le gradient en x � la valeur du pixel de l'image
			// autrement on fait la diff�rence entre le pixel et son voisin afin d'avoir des valeurs �lev�es aux contours dans l'image
			if (j + 1 >= sizeX) grady[j + sizeX*i] = imgIn[j + sizeX*i];
			else { grady[j + sizeX*i] = imgIn[j + 1 + sizeX*i] - imgIn[j + sizeX*i]; }
		}
	}

	// Afin d'obtenir le gradient de l'image, il reste plus qu'a faire la norme de l'addition des 2 matrices de gradient selon les diff�rents axes
	for (int i = 0; i<sizeY; i++) {
		for (int j = 0; j<sizeX; j++) {
			imgOut[j + sizeX*i] = sqrt(pow(gradx[j + sizeX*i], 2) + pow(grady[j + sizeX*i], 2));
		}
	}
}

// Fonction permettant de renvoyer le minimum entre 3 doubles
double min(double a, double b, double c){
	if(a<=b && a<=c)return a;
	else if (b<=c) return b;
	else return c;
}

// Fonction permettant de renvoyer le minimum entre 2 doubles
double min(double a, double b){
	return a <= b ? a : b;
}

// Recherche du chemin d'energie minimale dans imgIn 
void min_path(unsigned char *imgIn, int *posX, int *posY, int sizeX, int sizeY) 
{
	// Initialisation d'un tableau qui contiendra l'energie cumul�e de chaque pixel de l'energie de l'image et qui permettra de trouver le chemin le plus court pour la coupure
	int nrjcum[sizeX*sizeY];

	// Initialisation de la premi�re ligne de l'�nergie cumul�e � la premi�re ligne de valeurs de l'energie de l'image
	for(int j=0;j<sizeX;j++){
		nrjcum[j]=imgIn[j];		
	}			

	// Parcours sur toutes les lignes restantes � compl�ter de l'energie cumul�e
	for(int i=1; i<sizeY;i++){
		for(int j=0; j<sizeX;j++){
			
		// Si la colonne est la premiere de la ligne alors on cherche la valeur d'energie la plus faible entre les 2 valeurs d'�nergies en i-1,j et en i-1,j+1 de la ligne du dessus,
		// on l'additionne � la valeur d'energie � ce pixel de la matrice d'�nergie afin d'obtenir la valeur de l'energie cumul�e de ce pixel
		if(j-1<0){
			nrjcum[j+sizeX*i]= imgIn[j+sizeX*i]+min(nrjcum[j+sizeX*(i-1)],nrjcum[j+1+sizeX*(i-1)]);
				}
		// Si la colonne est la derniere de la ligne alors on cherche la valeur d'energie la plus faible entre les 2 valeurs d'�nergies en i-1,j-1 et en i-1,j de la ligne du dessus,
		// on l'additionne � la valeur d'energie � ce pixel de la matrice d'�nergie afin d'obtenir la valeur de l'energie cumul�e de ce pixel
		if(j+1>=sizeX){
				nrjcum[j+sizeX*i]= imgIn[j+sizeX*i]+min(nrjcum[j-1+sizeX*(i-1)],nrjcum[j+sizeX*(i-1)]);
				}
		// Si la colonne est situ� sur aucun des bords gauche ou droite alors on cherche la valeur d'energie la plus faible entre les 3 valeurs d'�nergies en i-1,j-1, en i-1,j et en i-1,j+1
		// de la ligne du dessus, on l'additionne � la valeur d'energie � ce pixel de la matrice d'�nergie afin d'obtenir la valeur de l'energie cumul�e de ce pixel
		else{
				nrjcum[j+sizeX*i]= imgIn[j+sizeX*i]+min(nrjcum[j-1+sizeX*(i-1)],nrjcum[j+sizeX*(i-1)],nrjcum[j+1+sizeX*(i-1)]);
			}
		}
	}

	//*********************************************************************
	// **Recherche de la derni�re valeur du chemin de plus petite �nergie**

	// Initialisation d'une variable avec la premi�re valeur de la derni�re ligne de notre �nergie cumul�e
	double temp=nrjcum[(sizeY-1)*sizeX];
	// Initialisation d'une variable contenant l'indice ou l'energie cumul�e de la ligne est minimum
	double mini = 0;

	// Parcours de la derniere ligne de l'energie cumul�e en commencant � l'indice 1 car nous avons deja recup�rer la valeur � l'indice 0 
	for(int j=1;j<sizeX;j++){
		// Si l'energie cumul�e du pixel est inferieur � celui enregistr� comme �tant le minimum de la ligne
		if(nrjcum[j+(sizeY-1)*sizeX]<temp) 
		{
			// On met � jour l'indice mini � j
			mini = j;
			// ainsi que sa valeur associ�e
			temp = nrjcum[j+(sizeY-1)*sizeX];
		}
	}
	//*********************************************************************

	double val;

	// Initialisation de la derniere valeur de chemin � mini puisque la derniere valeur minimum trouv�e est celle de la derni�re ligne
	posY[sizeY-1]= mini;

	// Parcours pour chaque ligne de l'energie cumul� (en partant de la derni�re ligne) des 3 valeurs voisines de la lignes de dessus 
	for(int i=sizeY-1; i>0;i--){
		
		// Si la valeur mini de la ligne actuelle est situ� sur le bord gauche alors on cherche la valeur d'�nergie la plus faible parmis les 2 energies cumul�es de la ligne de dessus
		// en i-1, pos et en i-1, pos+1 
		if(mini-1<0){
			val = min(nrjcum[posY[i]+sizeX*(i-1)],nrjcum[posY[i]+1+sizeX*(i-1)]);
			// Si la valeur mini est celle de droite alors on d�cale notre indice de mini de 1 vers la droite
			// sinon il reste au meme indice
			if(val == nrjcum[posY[i]+1+sizeX*(i-1)])
				mini = mini+1;
			}

		// Si la valeur mini de la ligne actuelle est situ� sur le bord droit alors on cherche la valeur d'�nergie la plus faible parmis les 2 energies cumul�es de la ligne de dessus
		// en i-1, pos-1 et en i-1, pos
		else if(mini+1>sizeX-1){
			val = min(nrjcum[posY[i]-1+sizeX*(i-1)],nrjcum[posY[i]+sizeX*(i-1)]);
			// Si la valeur mini est celle de gauche alors on d�cale notre indice de mini de 1 vers la gauche
			// sinon il reste au meme indice
			if(val == nrjcum[posY[i]-1+sizeX*(i-1)])
				mini = mini-1;
			}

		// Si la valeur mini de la ligne actuelle est situ� sur aucun des bords droit ou gauche alors on cherche la valeur d'�nergie la plus faible parmis les 3 energies cumul�es de la ligne de dessus
		// en i-1, pos-1, en i-1, pos et en i-1,pos+1
		else{
			val = min(nrjcum[posY[i]-1+sizeX*(i-1)],nrjcum[posY[i]+sizeX*(i-1)],nrjcum[posY[i]+1+sizeX*(i-1)]);
			// Si la valeur mini est celle de gauche alors on d�cale notre indice de mini de 1 vers la gauche
			// sinon il reste au meme indice
			if(val == nrjcum[posY[i]-1+sizeX*(i-1)])
				mini = mini-1;
			// Si la valeur mini est celle de droite alors on d�cale notre indice de mini de 1 vers la droite
			// sinon il reste au meme indice
			if(val == nrjcum[posY[i]+1+sizeX*(i-1)])
				mini = mini+1;
			}
		// On enregistre l'indice du mini � emprunter dans le chemin de la plus petite energie dans le tableau posY qui contiendra ce chemin
		posY[i-1]=mini;
		
	}
}

// Fonction permettant de retirer une coupure dans une image
void retirer_coupure(unsigned char *imgIn,int * posY, int sizeX, int sizeY, unsigned char *imgOut){
	
	// Parcours de chaque pixel de l'image d'entr�e
	for(int i=0; i<sizeY;i++){
		for(int j=0; j<sizeX;j++){
			// Lorsque nous avons d�pass� l'indice de la couture � enlever, nous �crasons les valeurs de la colonne precedante avec celle qu'il y a apres la cupure � enlever
			// ainsi de suite jusqu'a la fin de la ligne
			if(j>posY[i])
			{
				imgOut[(j-1)+i*(sizeX-1)]=imgIn[j+i*sizeX];
			}
			// Si nous ne sommes toujours pas arriv� � l'indice de la coupure sur la ligne alors on recopie simplement les pixels de l'image dans l'image de sortie
			else{
				imgOut[j +i*(sizeX-1)] =imgIn[j+i*sizeX];
			}
		}
	}
}

// Fonction permettant de recopier une tableau dans un autre
void recopie(unsigned char * imgIn, int sizeX, int sizeY, unsigned char * imgOut)
{
	for(int i=0;i<sizeY;i++)
	{
		for(int j=0;j<sizeX;j++)
		{
			imgOut[i*sizeX+j] = imgIn[i*sizeX+j];
		}
	}
}

// Fonction permettant de creer la transposer de l'image d'entr�e dans une l'image de sortie
void transpose(unsigned char * imgIn, int sizeX, int sizeY, unsigned char * imgOut)
{
	for(int i=0; i<sizeY;i++)
	{
		for(int j=0;j<sizeX;j++)
		{
			imgOut[j*sizeY+i] = imgIn[i*sizeX+j];
		}
	}
}

// Fonction permettant de creer l'histogramme d'une image
void histogramme(unsigned char *I, unsigned int* histo, int sizeX, int sizeY)
{
	for (int i=0; i<sizeY; i++){
		for(int j=0; j<sizeX; j++){
			histo[I[j+sizeX*i]]++;		
		}	
	}
}

// Fonction permettant de calculer la valeur  d'entropie d'un pixel en i,j de l'image d'entr�e avec son histogramme associ�
double entropy(unsigned char * imgIn, int i, int j,int sizeX,int sizeY,unsigned int * histo)
{
	unsigned int h1,w1;
	// Initialisation de la hauteur et de la largeur de l'�l�ment structurant
	h1=9;
    w1=9;
    	
    int N,M;

	// On initialise N et M � la moiti� de la hauteur et largeur de l'ES, ces valeurs nous permetteront de trouver
	// les coordonn�es des pixels de l'image
    N=h1/2;
    M=w1/2;
    	
    int x,y;
	double result = 0;
	// Pour chaque pixels de l'image on parcours chaque valeurs de l'ES
	for (int s=0; s < h1; s++){
    	for (int t = 0; t < w1; t++){
			// On r�cup�re les coordonn�es de pixels de l'image associ�es aux coordonn�es actuelles de l'ES(s,t)
			// pour r�cup�rer les valeurs de pixels de l'image point�es par ces coordonn�es.
    		x = i-N+s;
    		y = j-M+t;

			// Appel de la fontion g�rant les effets de miroir sur les bords de l'image
			miroir(x, y, imgIn,sizeX,sizeY);
    		
			// Ensuite on applique au r�sultat le calcul pour obtenir la valeur de l'entropie au pixel donn�
			float p_x= (float)histo[imgIn[y+sizeX*x]]/(sizeX*sizeY);
			if (p_x>0) result -= p_x*log(p_x)/log(2); 	
		}	
	}	
    return result;
}


void calcul_energie(unsigned char *imgIn, unsigned char *imgOut, int sizeX, int sizeY)
{
	
	// Initialisation d'un tableau qui contiendra l'histogramme de l'image
	unsigned int histogram[255];
	// Initialisation du tableau � 0
	for (int i=0; i<255;i++){
		histogram[i]=0;
	}
	// Remplissage de l'histogramme
	histogramme(imgIn,histogram,sizeX,sizeY);

	// Cr�ation d'une image contenant � chaque pixel la valeur de l'entropie du pixel associ� � l'image d'entr�e
	double entropie[sizeX*sizeY];
	// Parcours de chaque pixel de l'entropie/l'image d'entr�e
	for (int i=0; i < sizeY; i++){
    		for (int j = 0; j < sizeX; j++){
    		entropie[j+sizeX*i] = entropy(imgIn,i,j,sizeX,sizeY,histogram);
		}
	}

	// Calcul du gradient selon les 2 axes x et y de la m�me mani�re que dans la fonction IMgradient 
	double gradx[sizeX*sizeY];
	double grady[sizeY*sizeX];
	
	for(int i=0; i<sizeY; i++){
		for (int j=0; j<sizeX; j++){
			if(i+1>=sizeY) gradx[j+sizeX*i]=0;
			else {gradx[j+sizeX*i]=imgIn[j+sizeX*(i+1)]-imgIn[j+sizeX*i];}	
			
			if(j+1>=sizeX) grady[j+sizeX*i]=0;
			else {grady[j+sizeX*i]=imgIn[j+1+sizeX*i]-imgIn[j+sizeX*i];}
		}
	}

	// Enfin on applique la transformation du gradient sur chaque pixel avec chaque pixel de l'entropie afin d'obtenir un energie de l'image cr�er grace � une entropie
	for(int i =0; i<sizeY;i++){
		for(int j=0; j<sizeX;j++){
			imgOut[j+sizeX*i]=0.4*sqrt(pow(gradx[j+sizeX*i],2)+pow(grady[j+sizeX*i],2))+0.7*entropie[j+sizeX*i];
		}
	}

}


int main(int argc, char **argv)
{
 
	// Dimensions de l'image loup.pgm : 425,290
	// Dimensions de l'image mer.pgm : 353,500
	// Dimensions de l'image lena.pgm : 512,512
	int h =353;
	int w =500;
	unsigned char  * Iout = new unsigned char [h*w];

  // Lecture de l'image � redimensionner

 	char fileName[250];
	strcpy(fileName, "./images/mer.pgm");
	
	unsigned char  * Iori = new unsigned char [h*w];
	printf("\n Ouverture de %s de taille [%d,%d]", fileName, w, h);
	readPGM_Picture(fileName, Iori, w, h);

	//*******************************************
	//**********TRANSPOSE L'IMAGE**************** 
	unsigned char * It = new unsigned char[w*h];
	transpose(Iori,w,h,It);
	strcpy(fileName,"./images/transpose.pgm");
  	writePGM_Picture(fileName, It,h,w);
	int temp = h;	h=w;	w=temp;
	Iori = It; // On fait pointer l'image vers la transpos�e
	//*******************************************

	// Lecture d'un masque si besoin
	/*strcpy(fileName, "./images/mask1.pgm");

	unsigned char  * mask = new unsigned char [h*w];

	printf("\n Ouverture de %s de taille [%d,%d]", fileName, w, h);
	readPGM_Picture(fileName, mask, w, h);*/

	// Demande du type de calcul d'�nergie � effectuer
	int choix=0;
  	while(choix!=2&&choix!=1) {
		std::cout<<std::endl<<"Selection du type de calcul d'energie :" <<std::endl ;
		std::cout << "1. Gradient" <<std::endl ;
		std::cout << "2. Entropy"<<std::endl;
		std::cin >> choix ;
	}   
		
	// Demande du nombre de coupure selon les 2 directions
	int nb_coupure_vert, nb_coupure_horiz;
	std::cout<<std::endl<<"Nombre de coupure verticale : " <<std::endl ;
	std::cin>> nb_coupure_vert;
	std::cout<<std::endl<<"Nombre de coupure horizontale :" <<std::endl ;
	std::cin>> nb_coupure_horiz;

	int posX[w];
	int posY[h];
	
	unsigned char  * imgOut = new unsigned char [h*w];
	unsigned char  * imgOutMask = new unsigned char[h*w]; // Image permettant de stocker le resultat lorsque l'on retire une coupure sur le masque

  // Retirer coupure horizontale � plusieurs reprise
	for(int i = 0; i<nb_coupure_horiz;i++,w--)
	{		
	  // 1) Calculer la matrice d'�nergie de l'image 
		if(choix ==1)
			IMgradient(Iori,Iout,w,h);
			//IMgradient_mask(Iori,Iout,w,h,mask);
		else
			calcul_energie(Iori,Iout,w,h);
		
		printf("\n Ouverture de %s de taille [%d,%d]", fileName, w, h);
		strcpy(fileName,"./images/grad.pgm");
	  	writePGM_Picture(fileName, Iout,w,h);

	  // 2) D�terminer le chemin d'�nergie minimale 
		min_path(Iout, posX, posY, w, h);

	  // 3) L'enlever de l'image 
		retirer_coupure(Iori,posY, w, h, imgOut);
		strcpy(fileName,"./images/retire.pgm");
	  	writePGM_Picture(fileName, imgOut,w-1,h);

		//***********************************************
		//**********Utilisation d'un masque**************
		/*retirer_coupure(mask,posY, w, h, imgOutMask);
	
		strcpy(fileName,"./images/mask_modif.pgm");
	  	writePGM_Picture(fileName, imgOutMask,w-1,h);

		readPGM_Picture(fileName, mask, w-1, h);*/
		//***********************************************

		// On recopie l'image obtenu dans l'image d'entr�e (qui va subir les modifications)
		recopie(imgOut, w, h, Iori);
	}
	
	// Dans le cas ou on a choisi de ne pas faire de coupure horizontale on recopie l'image dans l'image de sortie afin d'effectuer les modifications dessus
	if(nb_coupure_horiz == 0)
		recopie(Iori,w,h,imgOut);

	//*******************************************
	//**********TRANSPOSE L'IMAGE**************** 
	//It = new unsigned char[w*h];
	transpose(imgOut,w,h,It);
	strcpy(fileName,"./images/retire.pgm");
  	writePGM_Picture(fileName, It,h,w);
	temp = h;	h=w;	w=temp;
	Iori = It; // On fait pointer l'image vers la transpos�e
	//*******************************************

	// Retirer coupure verticale � plusieurs reprise
	for(int i = 0; i<nb_coupure_vert;i++,w--)
	{		

	  // 1) Calculer la matrice d'�nergie de l'image 
		if(choix ==1)
			IMgradient(Iori,Iout,w,h);
			//IMgradient_mask(Iori,Iout,w,h,mask);
		else
			calcul_energie(Iori,Iout,w,h);

		printf("\n Ouverture de %s de taille [%d,%d]", fileName, w, h);
		strcpy(fileName,"./images/grad.pgm");
	  	writePGM_Picture(fileName, Iout,w,h);

	  // 2) D�terminer le chemin d'�nergie minimale 
		min_path(Iout, posX, posY, w, h);

	  // 3) L'enlever de l'image 
		retirer_coupure(Iori,posY, w, h, imgOut);
		strcpy(fileName,"./images/retire.pgm");
	  	writePGM_Picture(fileName, imgOut,w-1,h);

		//***********************************************
		//**********Utilisation d'un masque**************
		/*retirer_coupure(mask,posY, w, h, imgOutMask);
	
		strcpy(fileName,"./images/mask_modif.pgm");
	  	writePGM_Picture(fileName, imgOutMask,w-1,h);

		readPGM_Picture(fileName, mask, w-1, h);*/
		//***********************************************

		// On recopie l'image obtenu dans l'image d'entr�e (qui va subir les modifications)
		recopie(imgOut, w, h, Iori);
	}

    return 0;
}
