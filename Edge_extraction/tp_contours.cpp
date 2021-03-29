/****************************************************************************
 * Copyright (C) 2016 Universite de Rennes 1. All rights reserved.
 *
 * This software was developed at:
 * Universite de Rennes 1
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 *
 * This file uses the ViSP library.
 *****************************************************************************/

/****************************************************************************
 * NOMS - PRENOMS:
 *  - GLEDEL Clément
 *	- DUGNE Quentin
 * 
 * Date :
 *****************************************************************************/


#include <iostream>

#include <visp/vpConfig.h>
#include <visp/vpDebug.h>

#include <visp/vpImage.h>
#include <visp/vpImageFilter.h>
#include <visp/vpImageIo.h>
#include <visp/vpDisplayX.h>

using namespace std;


void afficheImage(vpImage<unsigned char> img, int posX, int posY, const char *title)
{
    vpDisplayX d(img, posX, posY, title);
    vpDisplay::display(img);
    vpDisplay::flush(img);
    vpDisplay::getClick(img);
    vpDisplay::close(img);
}

/**
 * @brief calcule la rotation de k x pi/4 du gabarit
 * @param gabarit: gabarit d'entrée
 * @param k: parametre de la rotation (k x pi/4)
 * @return gabarit apres rotation de k x pi/4 (nouvelle image, allocation memoire)
 */
vpImage<int> rotation_masque(const vpImage<int> &gabarit, int k)
{
    
    vpImage<int> rotated(gabarit.getHeight(),gabarit.getWidth(),0) ;
    
    int half_h=gabarit.getHeight()/2;
    int half_w=gabarit.getWidth()/2;
    int offseti, offsetj;
    
    for(int i=-half_h;i<=half_h;i++){
        for(int j=-half_h;j<=half_h;j++){
            offseti=round(i*cos(k*M_PI/4)-j*sin(k*M_PI/4));
            offsetj=round(i*sin(k*M_PI/4)+j*cos(k*M_PI/4));
            rotated[offseti+half_h][offsetj+half_h]=gabarit[i+half_h][j+half_h];
        }
    }
    return rotated;
}

// Fonction permettant de gérer les effets de miroir sur les bords de l'image
void miroir(int &x, int &y, const vpImage<unsigned char > &I) {

	// Si la valeur de x est inférieur à 0 et est donc en dehors (au dessus visuellement) de l'image on prend 
	// la valeur récupérée par effet de miroir par rapport au bord du dessus de l'image
	if (x<0) {
		x = -x - 1; // on met la valeur négatif -> positive puis on soustrait 1 car la valeur 0 fait partie de l'image (bord)					
	}
	// Si la valeur de y est inférieur à 0 et est donc en dehors (à gauche visuellement) de l'image on prend 
	// la valeur récupérée par effet de miroir par rapport au bord gauche de l'image
	if (y<0) {
		y = -y - 1; // on met la valeur négatif -> positive puis on soustrait 1 car la valeur 0 fait partie de l'image (bord)
	}
	// Si la valeur de x est supérieur ou égale à la hauteur de l'image et est donc en dehors (en dessous visuellement)
	// de l'image on prend la valeur récupérée par effet de miroir par rapport au bord du dessous de l'image
	if (x >= I.getHeight()) {
		x = I.getHeight() - 1 - (x - I.getHeight()); // La valeur au bord moins la différence entre le bord et la valeur de x
	}
	// Si la valeur de y est supérieur ou égale à la largeur de l'image et est donc en dehors (à droite visuellement)
	// de l'image on prend la valeur récupérée par effet de miroir par rapport au bord droite de l'image
	if (y >= I.getWidth()) {
		y = I.getWidth() - 1 - (y - I.getWidth()); // La valeur au bord moins la différence entre le bord et la valeur de y
	}
}

// Fonction permettant de gérer les effets de miroir sur les bords de l'image
void miroir(int &x, int &y, const vpImage<double > &I) {

	// Si la valeur de x est inférieur à 0 et est donc en dehors (au dessus visuellement) de l'image on prend 
	// la valeur récupérée par effet de miroir par rapport au bord du dessus de l'image
	if (x<0) {
		x = -x - 1; // on met la valeur négatif -> positive puis on soustrait 1 car la valeur 0 fait partie de l'image (bord)					
	}
	// Si la valeur de y est inférieur à 0 et est donc en dehors (à gauche visuellement) de l'image on prend 
	// la valeur récupérée par effet de miroir par rapport au bord gauche de l'image
	if (y<0) {
		y = -y - 1; // on met la valeur négatif -> positive puis on soustrait 1 car la valeur 0 fait partie de l'image (bord)
	}
	// Si la valeur de x est supérieur ou égale à la hauteur de l'image et est donc en dehors (en dessous visuellement)
	// de l'image on prend la valeur récupérée par effet de miroir par rapport au bord du dessous de l'image
	if (x >= I.getHeight()) {
		x = I.getHeight() - 1 - (x - I.getHeight()); // La valeur au bord moins la différence entre le bord et la valeur de x
	}
	// Si la valeur de y est supérieur ou égale à la largeur de l'image et est donc en dehors (à droite visuellement)
	// de l'image on prend la valeur récupérée par effet de miroir par rapport au bord droite de l'image
	if (y >= I.getWidth()) {
		y = I.getWidth() - 1 - (y - I.getWidth()); // La valeur au bord moins la différence entre le bord et la valeur de y
	}
}

// Fonction permettant d'effectuer la convolution d'une image I par une matrice ligne ou colonne Kx
void filtrage1D(const vpImage< double > &I, vpImage< double > &Ic, const vpMatrix &Kx)
{
	 unsigned int  h,w;
   
    h=I.getHeight(); //hauteur de l'image 
    w=I.getWidth(); //largeur de l'image
    
    unsigned int h1,w1;
    
    h1= Kx.getRows(); //hauteur du masque
    w1= Kx.getCols(); //largeur du masque
    
    unsigned int dim = h1; // dim correspond à la valeur la plus grande entre largeur et hauteur du masque
    if (h1<w1) dim=w1; //
    
    int N,M; 
    	
    N=h1/2; //N égal à la hauteur du filtre divisé par 2
    M=w1/2; //M égal à la largeur du filtre divisé par 2
    	
    int x,y; //coordonnées utile pour faire le miroir
 
    //Allocate memory for an [height x width] image and initialize the image to val.
    Ic.resize(I.getHeight(),I.getWidth(),0) ;
    
	// Parcours sur chacun des pixels de l'image
    for (int i=0; i<h; i++){
    	for(int j =0; j<w; j++){
			// Parcours sur toutes les valeurs de la matrice de masquage
    		for(int k=0; k<dim; k++){   		
				// Si la matrice est une matrice ligne
    			if(h1<w1){
					// Selon le même principe que la fonction miroir ici on traite seulement le cas où la valeur de coordonnée
					// sort de l'image selon l'axe des y (colonne) car c'est une matrice ligne
					y = j - M + k;
					if (y<0) {
						y = -y - 1;
					}
					if (y >= I.getWidth()) {
						y = I.getWidth() - 1 - y + I.getWidth();
					}
					// On multiplie la valeur de pixel de l'image à la valeur de la matrice associée que l'on ajoute à la valeur du pixel de la nouvelle image résultante
					Ic[i][j] += Kx[0][k]*I[i][y];
    			}
    			// Si la matrice est une matrice colonne
    			else if(h1>w1){
					// Selon le même principe que la fonction miroir ici on traite seulement le cas où la valeur de coordonnée
					// sort de l'image selon l'axe des x (ligne) car c'est une matrice ligne
					x = i - N + k;
					if (x<0) {
						x = -x - 1;
					}
					if (x >= I.getHeight()){
						x = I.getHeight() - 1 - (x - I.getHeight());
					}
					// On multiplie la valeur de pixel de l'image à la valeur de la matrice associée que l'on ajoute à la valeur du pixel de la nouvelle image résultante
    				Ic[i][j] += Kx[k][0]*I[x][j];
    			}
    		}
    	}
    }
}

/**
* Filtrage separable : K= Kx * Ky => I*K = (I*Kx)*Ky
 * avec :
 * K de taille N x N
 * Kx de taille 1 x N
 * Ky de taille N x 1
 * I*K = t[t(I*Kx)*tKy] avec t= transposition, Kx, Ky de taille 1 x N
 **/
void filtrage2D_separable(const vpImage< unsigned char > &I, vpImage< double > &Ic, const vpMatrix &Kx, const vpMatrix &Ky)
{
	// Image intermédiaire utilisé en tant que conversion de l'image d'entrée en vpImage<double>
    vpImage< double > Id;
    vpImageConvert::convert(I,Id); // convert unsigned char to double (with memory allocation)

	unsigned int h, w;
	h=I.getHeight(); //hauteur de l'image 
    w=I.getWidth(); //largeur de l'image

	// Image intermédiaire utilisé pour récupérer le résultat de la convolution selon la matrice Kx
	vpImage<double> Iinter(h,w);
	// Premier appel du filtrage selon une matrice ligne ou colone Kx
	filtrage1D(Id, Iinter, Kx);
	// Second appel du filtrage selon une matrice ligne ou colone Ky
	filtrage1D(Iinter, Ic, Ky);
    
    Id.destroy();
}



// Fonction permettant d'effectuer la convolution d'un pixel d'une image I au coordonées [i,j] par une matrice K
double convolution(const vpImage< unsigned char > &I, int i, int j, const vpImage<int> &K)
{
	unsigned int h1,w1;
	// Récupération de la hauteur et de la largeur de la matrice de convolution
	h1=K.getRows();
    w1=K.getCols();
    	
    int N,M;
    // On initialise N et M à la moitié de la hauteur et largeur de la matrice, ces valeurs nous permetteront de trouver
	// les coordonnées des pixels de l'image
    N=h1/2;
    M=w1/2;
    	
    int x,y;
    double result = 0;
    
	// Pour chaque pixels de l'image on parcours chaque valeurs de la matrice
    for (int s=0; s < h1; s++){
    	for (int t = 0; t < w1; t++){
			// On récupère les coordonnées de pixels de l'image associées aux coordonnées actuelles de la matrice
			// pour récupérer les valeurs de pixels de l'image pointées par ces coordonnées.
    		x = i-N+s;
    		y = j-M+t;

			// Appel de la fontion gérant les effets de miroir sur les bords de l'image
			miroir(x, y, I);
    		
			// On multiplie la valeur de pixel de l'image à la valeur de la matrice associée que l'on ajoute au résultat
    		result += K[s][t]*I[x][y];	
		}
	}  
    return result;
}

// Fonction permettant de retourner le gradient de l'image im en utilisant la méthode de Kirsh
void gradient(const vpImage< unsigned char > & im, vpImage< double > & Result)
{
	// Récupération de la hauteur et largeur de l'image 
	int h = im.getHeight();
	int w = im.getWidth();

	Result.resize(h, w, 0);
	
	// Création d'une image qui contiendra les max du gradient selon toutes les directions
	vpImage <int> Gradmax(h, w, 0);
	// Création d'une image qui contiendra les directions dans lesquelles les max ont été trouvés pour chaque pixel de l'image
	vpImage <int> Orientation(h, w, 0);
	// Gabarit détenant la première orientation du masque
	vpImage <int> gabarit(3, 3, 0);
	gabarit[0][0] = -3; gabarit[0][1] = -3; gabarit[0][2] = 5;
	gabarit[1][0] = -3; gabarit[1][1] = 0; gabarit[1][2] = 5;
	gabarit[2][0] = -3; gabarit[2][1] = -3; gabarit[2][2] = 5;

	// Images intermediaires
	vpImage <int> gab_intermediaire(3, 3, 0);
	vpImage <int> inter(h,w,0);
		
	// Parcours pour chaque rotation du masque
    for (int k=0; k<8; k++){
		// Recupération de la rotation du masque d'un angle de k*pi/4
		gab_intermediaire=rotation_masque(gabarit,k);
		
		// Parcours de chaque pixel de l'image
		for (int i=0; i<h; i++){
			for(int j=0; j<w; j++){
				// Remplissage de l'image intermédiaire par la convolution de l'image d'entrée avec la rotation du masque d'un angle de k*pi/4
				inter[i][j] = convolution(im,i,j, gab_intermediaire);
				// Si la valeur renvoyé par la convolution de l'image est supérieur à la valeur présente dans l'image des max
				if(abs(inter[i][j])>Gradmax[i][j])
				{
					// Enregistrement de la nouvelle valeur de max du pixel dans l'image des max
					Gradmax[i][j] = inter[i][j];
					// Enregistrement de la nouvelle direction détenant ce max
					Orientation[i][j]= k;
				}	
			}
		}
	}

	// Parcours de chaque pixel de l'image (sans les bordures)
	for (int i=1; i<h-1; i++){
		for(int j=1; j<w-1; j++){
			// Récupération de l'orientation au pixel
			int k = Orientation[i][j];

			// Si cette orientation est orienté de gauche à droite ou inversement
			if(k==0 || k ==4){
				// Si une des 2 valeurs voisines selon cette orientation est supérieur alors ce n'est pas un contour (on met la valeur à 0)
				if(Gradmax[i][j-1]>=Gradmax[i][j] || Gradmax[i][j+1]>=Gradmax[i][j])
				{
					Result[i][j] = 0;
				}
				// Sinon (elle est supérieur au 2 voisines) c'est donc un contour
				else Result[i][j]=Gradmax[i][j];
			}

			// Si cette orientation est orienté de haut/droite à bas/gauche ou inversement
			else if(k==1 || k ==5){
				// Si une des 2 valeurs voisines selon cette orientation est supérieur alors ce n'est pas un contour (on met la valeur à 0)
				if(Gradmax[i-1][j+1]>=Gradmax[i][j] || Gradmax[i+1][j-1]>=Gradmax[i][j])
				{
					Result[i][j] = 0;
				}
				// Sinon (elle est supérieur au 2 voisines) c'est donc un contour
				else Result[i][j]=Gradmax[i][j];	
			}

			// Si cette orientation est orienté de haut en bas ou inversement
			else if(k==2 || k ==6){
				// Si une des 2 valeurs voisines selon cette orientation est supérieur alors ce n'est pas un contour (on met la valeur à 0)
				if(Gradmax[i-1][j]>=Gradmax[i][j] || Gradmax[i+1][j]>=Gradmax[i][j])
				{
					Result[i][j] = 0;
				}
				// Sinon (elle est supérieur au 2 voisines) c'est donc un contour
				else Result[i][j]=Gradmax[i][j];
			}

			// Si cette orientation est orienté de haut/gauche à bas/droite ou inversement
			else if(k==3 || k ==7){
				// Si une des 2 valeurs voisines selon cette orientation est supérieur alors ce n'est pas un contour (on met la valeur à 0)
				if(Gradmax[i-1][j-1]>=Gradmax[i][j] || Gradmax[i+1][j+1]>=Gradmax[i][j])
				{
					Result[i][j] = 0;
				}
				// Sinon (elle est supérieur au 2 voisines) c'est donc un contour
				else Result[i][j]=Gradmax[i][j];
			}
		}
	}			
}

// Fonction permettant d'effectuer un Zero-Crossing d'une image
void zero_crossing(const vpImage< double > & inter, vpImage< double > & Result, double seuil)
{
	// Récupération de la hauteur et largeur de l'image 
	int h = inter.getHeight();
	int w = inter.getWidth();

	double max = 0;
	// Parcours de chaque pixel de l'image
	for (int i = 1; i<h - 1; i++) {
		for (int j = 1; j<w - 1; j++) {

			// Si les voisins haut/gauche et bas/droite sont négatif et positif ou inversement
			if ((inter[i - 1][j - 1]<0 && inter[i + 1][j + 1]>0) || (inter[i - 1][j - 1]>0 && inter[i + 1][j + 1]<0)) {
				// Test si la différence est supérieur au max enregistré si oui alors on prend cette nouvelle valeur comme max
				if (abs(inter[i - 1][j - 1] - inter[i + 1][j + 1])>max) 
					max = abs(inter[i - 1][j - 1] - inter[i + 1][j + 1]);
			}

			// Si les voisins haut et bas sont négatif et positif ou inversement
			if ((inter[i - 1][j]<0 && inter[i + 1][j]>0) || (inter[i - 1][j]>0 && inter[i + 1][j]<0)) {
				// Test si la différence est supérieur au max enregistré si oui alors on prend cette nouvelle valeur comme max
				if (abs(inter[i - 1][j] - inter[i + 1][j])>max) 
					max = abs(inter[i - 1][j] - inter[i + 1][j]);
			}

			// Si les voisins gauche et droite sont négatif et positif ou inversement
			if ((inter[i][j - 1]<0 && inter[i][j + 1]>0) || (inter[i][j - 1]>0 && inter[i][j + 1]<0)) {
				// Test si la différence est supérieur au max enregistré si oui alors on prend cette nouvelle valeur comme max
				if (abs(inter[i][j - 1] - inter[i][j + 1])>max) 
					max = abs(inter[i][j - 1] - inter[i][j + 1]);
			}

			// Si les voisins haut/droite et bas/gauche sont négatif et positif ou inversement
			if ((inter[i - 1][j + 1]<0 && inter[i + 1][j - 1]>0) || (inter[i - 1][j + 1]>0 && inter[i + 1][j - 1]<0)) {
				// Test si la différence est supérieur au max enregistré si oui alors on prend cette nouvelle valeur comme max
				if (abs(inter[i - 1][j + 1] - inter[i + 1][j - 1])>max) 
					max = abs(inter[i - 1][j + 1] - inter[i + 1][j - 1]);
			}

			// Si la valeur max est supérieur au seuil alors c'est un contour
			if (max > seuil) Result[i][j] = 255;

			// Réinitialisation de max à 0
			max = 0;
		}
	}
}

// Fonction permettant d'obtenir une image de contour selon la méthode du Laplacien (méthode utilisant le filtrage 2D séparable vu au tp4 de BINP) et du zero crossing
void Laplacien(const vpImage< unsigned char > & im, vpImage< double > & Result, int pourcentage_seuil)
{
	// Récupération de la hauteur et largeur de l'image 
	int h = im.getHeight();
	int w = im.getWidth();

	Result.resize(h,w,0);
	
	// Image de stockage des image filtrées
	vpImage<double> Ifiltre1(h, w, 0);
	vpImage<double> Ifiltre2(h, w, 0);
	// Image intermédiaire
	vpImage<double> inter(h, w, 0);

	// Matrices ligne et colonne contenant les masque 1D de filtrage
	vpMatrix Kx(1, 3);
	vpMatrix Ky(3, 1);

	// 1er filtre de Laplacien
	Kx[0][0]=1;Kx[0][1]=1;Kx[0][2]=1;
	Ky[0][0]=-1;
	Ky[1][0]=2;
	Ky[2][0]=-1;

	// Application du filtre 2D séparable sur notre image avec la matrice ligne Kx et la matrice colonne Ky
	filtrage2D_separable(im, Ifiltre1, Kx, Ky);
	// Parcours de chaque pixel de l'image
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
		{
			// Normalisation des valeurs (division par 8 car on ajout la contribution de 8 voisins)
			Ifiltre1[i][j] = Ifiltre1[i][j]/8;
		}
	}

	// 2eme filtre de Laplacien
	Kx[0][0]=-1;Kx[0][1]=2;Kx[0][2]=-1;
	Ky[0][0]=1;
	Ky[1][0]=1;
	Ky[2][0]=1;

	// Application d'un deuxième filtrage 2D séparable sur notre image avec la nouvelle matrice ligne Kx et colonne Ky
	filtrage2D_separable(im, Ifiltre2, Kx, Ky);
	
	// Maximum de l'image résultante
	double maximum = 0;
	// Parcours de chaque pixel de l'image
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
		{
			// Normalisation des valeurs du deusième filtrage (division par 8 car on ajout la contribution de 8 voisins)
			// Puis ajout des valeurs du premier filtrage
			inter[i][j] = Ifiltre2[i][j]/8 + Ifiltre1[i][j];
			// Récupération de la valeur maximal dans l'image résultante
			if(inter[i][j] > maximum) maximum = inter[i][j];
		}
	}
	
	// Calcul du seuil en fonction du pourcentage donné en paramètre
	double seuil = maximum * pourcentage_seuil/100;
	
	// Appel de la fonction effectuant un zero crossing de l'image
	zero_crossing(inter, Result, seuil);
}

// Fonction permettant d'obtenir une image de contour en utilisant la méthode de LoG approximé par une Différence de Gaussiennes (DoG) et du zero crossing
// sigma1 est le paramètre de la première gausienne et sigma2 et le paramètre de la deuxième gaussienne
// le rapport sigma1/sigma2 devant être égale à 1,6 et dont les valeurs 6*sigma1+1 et 6*sigma2+1 doivent être impaires
void DOG(const vpImage< unsigned char > & im, vpImage< double > & Result, double sigma1, double sigma2, int pourcentage_seuil)
{
	// Récupération de la hauteur et de la largeur de l'image
	int h = im.getHeight();
	int w = im.getWidth();

	// Détermination des tailles des gausiennes voisines avec les sigmas donnés en paramètre
	int size1 = round(sigma1*6+1);
	int size2 = round(sigma2*6+1);

	// Filtre 1 qui contiendra la première gaussienne décomposée en une matrice ligne et une matrice colonne
	vpMatrix filtre1c(size1,1);
	vpMatrix filtre1l(1,size1);
	// Filtre 2 qui contiendra la deuxième gaussienne décomposée en une matrice ligne et une matrice colonne
	vpMatrix filtre2c(size2,1);
	vpMatrix filtre2l(1,size2);
	
	// Création des 2 tableaux qui contiendront les valeurs d'une partie des gaussiennes à symétriser
	double Gauss1[(size1+1)/2];
	double Gauss2[(size2+1)/2];
	
	// Images intermédiares
	vpImage <double> inter(h,w,0);
	vpImage <double> inter1(h,w,0);
	vpImage <double> inter2(h,w,0);
	
	// Récupération des valeurs des gaussiennes
	vpImageFilter::getGaussianKernel((double *) &Gauss1,size1);
	vpImageFilter::getGaussianKernel((double *) &Gauss2,size2);
	
	// Parcours de chaque valeur des tableaux contenant les valeurs de la première gaussienne
	for(int i=0; i<(size1)/2 +1;i++)
	{
		// Remplissage de la matrice colonne du filtre1 en prenant en compte la symétrie à faire
		filtre1c[i][0] = Gauss1[(size1)/2-i];
		filtre1c[size1-1-i][0] = Gauss1[(size1)/2-i];

		// Remplissage de la matrice ligne du filtre1 en prenant en compte la symétrie à faire
		filtre1l[0][i] = Gauss1[(size1)/2-i];
		filtre1l[0][size1-1-i] = Gauss1[(size1)/2-i];
	}

	// Parcours de chaque valeur des tableaux contenant les valeurs de la deuxième gaussienne
	for(int i=0; i<(size2)/2 +1;i++)
	{
		// Remplissage de la matrice colonne du filtre2 en prenant en compte la symétrie à faire
		filtre2c[i][0] = Gauss2[(size2)/2-i];
		filtre2c[size2-1-i][0] = Gauss2[(size2)/2-i];

		// Remplissage de la matrice ligne du filtre2 en prenant en compte la symétrie à faire
		filtre2l[0][i] = Gauss2[(size2)/2-i];
		filtre2l[0][size2-1-i] = Gauss2[(size2)/2-i];
	}
	// Récupération du filtrage de l'image par le filtre 1 dans l'image inter1
	filtrage2D_separable(im, inter1, filtre1c, filtre1l);
	// Récupération du filtrage de l'image par le filtre 2 dans l'image inter2
	filtrage2D_separable(im, inter2, filtre2c, filtre2l);
	
	double maximum = 0;
	// Parcours de chaque pixel de l'image
	for (int i = 0; i<h; i++){
		for(int j=0; j<w; j++)
		{
			// Remplissage de l'image filtré en faisant la soustraction des 2 images filtrées (principe du DoG)
			inter[i][j]= inter1[i][j]-inter2[i][j];
			// Récupération du max de l'image résultante
			if(inter[i][j] > maximum) maximum = inter[i][j];
		}
	}
	
	// Calcul du seuil en fonction du pourcentage donné en paramètre
	double seuil = maximum * pourcentage_seuil/100;

	// Appel de la fonction effectuant un zero crossing de l'image
	zero_crossing(inter, Result, seuil);
}

// Fonction permettant de savoir si le pixel de l'image im (détecté comme un contour) correspond à un contour dans l'image de référence dans un voisinage de 3x3
bool contour_correct(const vpImage< unsigned char > & im, const vpImage< unsigned char > & ref, int i, int j, int valeur_contour)
{
    int x,y;
    
	// Pour chaque pixels de l'image on parcours chaque valeurs de la matrice
    for (int s=0; s < 3; s++){
    	for (int t = 0; t < 3; t++){
			// On récupère les coordonnées de pixels de l'image associées aux coordonnées actuelles de la matrice (voisinage 3x3 dans notre cas)
			// pour récupérer les valeurs de pixels de l'image pointées par ces coordonnées.
    		x = i-1+s;
    		y = j-1+t;

			// Appel de la fontion gérant les effets de miroir sur les bords de l'image
			miroir(x, y, im);
    		
			// Si un contour est présent dans le voisinage 3x3 de l'image ref alors le contour est correct
			if(ref[x][y] == valeur_contour) return true;
		}
	}  
	return false;
}

void nb_contours(const vpImage< unsigned char > & im, const vpImage< unsigned char > & ref)
{
	// Initialisation des données que l'on cherche à obtenir
	int contour_det =0;
	int contour_ref =0;
	int corrects =0;
	int faux_positif =0;
	int faux_negatif =0;
	
	bool correct = false;
	
	// Parcours de chaque pixel de l'image
	for(int i=0; i<im.getHeight(); i++)
	{
		for(int j=0; j<im.getWidth(); j++)
		{
			// Si la valeur est un contour dans l'image de contour à tester
			if(im[i][j] ==255)
			{
				// Incrémente le nombre de contour détecté dans l'image à tester
				 contour_det++;
				 // On récupère le résultat du test de la présence d'un contour dans le voisinage 3x3 de l'image ref
				 correct = contour_correct(im,ref,i,j,0);
				 // S'il est correct alors on incrémente le nombre de contour correct
				 if(correct) corrects++;
			}

			// Si la valeur est un contour dans l'image de référence
			if(ref[i][j] ==0) 
			{
				// On incrémente le nombre de contour dans l'image de référence
				contour_ref++;
			}
		}
	}
	// Calcul du nombre de faux positif et de faux négatif
	faux_positif = contour_det - corrects;
	faux_negatif = contour_ref - corrects;
	
	// Calculs des valeurs de performance et de taux
	double performance = (double) corrects / (corrects + faux_positif + faux_negatif);
	double tauxfp = (double) faux_positif / (corrects + faux_positif + faux_negatif);
	double tauxfn = (double) faux_negatif / (corrects + faux_positif + faux_negatif);

	// Affichage des différentes valeurs
	cout << " Nb de contours : " << contour_det << endl ;
	cout << " Nb de contours ref : " << contour_ref << endl ;
	cout << " Nb de contours corrects : " << corrects << endl ;
	cout << " Nb de contours faux positif : " << faux_positif << endl ;
	cout << " Nb de contours faux négatif : " << faux_negatif << endl ;
	cout << " Performance : " << performance << endl ;
	cout << " Taux de faux positif : " << tauxfp << endl ;
	cout << " Taux de faux négatif : " << tauxfn << endl ;
}


int main(int argc, char **argv)
{
	
	cout << "TIA TP : DETECTION DE CONTOURS " << endl ;
	cout << "--" << endl ;
	
	// Image originale
	vpImage<unsigned char>  I0;
	vpImageIo::read(I0,"../../images/contours/gnu.pgm");//lions.pgm");//bear_2.pgm");
	afficheImage(I0,100,100,"Image originale");
  
	// Image contour de référence
	vpImage<unsigned char>  Icontour;
	vpImageIo::read(Icontour,"../../images/contours/gt/gnu_gt_binary.pgm");//lions_gt_binary.pgm");//bear_2_gt_binary.pgm");
	afficheImage(Icontour,100,100,"Image contour");
	
	// Initialisation mages qui contiendront les résulats
	vpImage<double> Igrad(I0.getHeight(),I0.getWidth(),0);
	vpImage<unsigned char> result(I0.getHeight(),I0.getWidth(),0);
	vpImage<double> ILaplacien(I0.getHeight(),I0.getWidth(),0);
	vpImage<unsigned char> result2(I0.getHeight(),I0.getWidth(),0);
	vpImage<double> IDOG(I0.getHeight(),I0.getWidth(),0);
	vpImage<unsigned char> result3(I0.getHeight(),I0.getWidth(),0);
	vpImage<unsigned char> seuillage(I0.getHeight(),I0.getWidth(),0);
	
	// Récupérations des résultats
	gradient(I0, Igrad);
	Laplacien(I0, ILaplacien, 50);
	DOG(I0, IDOG, 6.4,4,21);//3.68, 2.3, 38);
	
	// Convertion des résultats en image d'unsigned char
	vpImageConvert::convert(Igrad, result);
	vpImageConvert::convert(ILaplacien, result2);
	vpImageConvert::convert(IDOG, result3);
		
	// *****SEUILLAGE DU GRADIENT*****
	for (int i=0; i<I0.getHeight(); i++){
		for(int j=0; j<I0.getWidth(); j++){
		
			if(result[i][j]>89){
				seuillage[i][j]=255;
			}
			else seuillage[i][j]=0;
		}
	}
	// *******************************
	
	vpImageIo::write(result,"../../images/result_image_grad.pgm");
	afficheImage(result,100,100,"Image Gradient");
	
	cout << endl << " Gradient seuillé : " << endl << endl ;
	nb_contours(seuillage,Icontour);
	vpImageIo::write(seuillage,"../../images/result_image_grad_seuil.pgm");
	afficheImage(seuillage,100,100,"Image Gradient seuillé");
	
	cout << endl << " Laplacien : " << endl << endl ;
	nb_contours(result2,Icontour);
	vpImageIo::write(result2,"../../images/result_image_laplacien.pgm");
	afficheImage(result2,100,100,"Image laplacien");
	
	cout << endl << " DOG : " << endl << endl ;
	nb_contours(result3,Icontour);
	vpImageIo::write(result3,"../../images/result_image_DOG.pgm");
	afficheImage(result3,100,100,"Image DOG");
	
	cout << "Fin du programme " << endl ;

	I0.destroy();
	Icontour.destroy();
	Igrad.destroy();
	ILaplacien.destroy();
	IDOG.destroy();
	result.destroy();
	result2.destroy();
	result3.destroy();
	return(0);
}















