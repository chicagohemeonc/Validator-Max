#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include "peptide_fragment.hpp"
using namespace std;


PeptideFragment::PeptideFragment ( int massType ) {
  _iMassType = massType;
  _initialize_mass();
  _init_comp();
}


void PeptideFragment::analyze (char *peptide, int charge, char *modification) {
  unsigned int i, j;
  double dProton = PROTON_MASS;

  a.clear();
  b.clear();
  c.clear();
  x.clear();
  y.clear();
  z.clear();
  zdot.clear();
  mass.clear();

  /* remove non-alpha chars */
  j = 0;
  for (i = 0; i < strlen(peptide); i++)
    if (isalpha(peptide[i]))
      _szInputSequence[j++] = peptide[i];
  _szInputSequence[j] = '\0';
  _iLenPeptide = j;

  _iCharge = charge;
  _modification = string(modification);


  for (int i = 0; i < MAX_SEQUENCE; i++)
    pdPositionMod[i] = 0;

  if ( _modification.length() ) {
    if ( _modification == "silac" or _modification == "SILAC" ) {
      pdMassAA['K'] += 8.014199;
      pdMassAA['R'] += 10.008269;
    }
    else if ( _modification == "18O" ) {
      if (_iMassType==1) /*mono*/
        pdMassAA['R'] += 3.98814;
      else
        pdMassAA['R'] += 3.9737;
    }
    else if ( _modification == "CAM-C" or _modification == "carbamidomethyl C") {
      if (_iMassType==1) /*mono*/
        pdMassAA['C'] += 57.021464;
      else
        pdMassAA['C'] += 57.0513;
    }
    else {
      char *tok;
      tok = strtok( (char *)_modification.c_str(), " ");
      while (tok != NULL) {
        char *pStr;
        double dMass;
        char szResidue[100];

        pStr=strchr(tok, '@');
        *pStr = ' ';

        sscanf(tok, "%lf %s", &dMass, szResidue);

        /*
         * first check if residue or number is entered
         */
        if (strspn(szResidue, "0123456789")==strlen(szResidue)) {
          int iPos=0;

          sscanf(szResidue, "%d", &iPos);
          if (iPos>=1)
            pdPositionMod[iPos-1] = dMass;
        }
        else {
          pdMassAA[toupper(szResidue[0])] += dMass;
        }

        tok = strtok(NULL, " ");
      }
    }
  }

  dPepMass = pdMassAA['o'] + pdMassAA['h'] + pdMassAA['h'] + dProton;
  for ( i = 0; i < _iLenPeptide; i++ ) {
    dPepMass += pdMassAA[(int)(_szInputSequence[i])];
    dPepMass += pdPositionMod[i];
  }}



std::string PeptideFragment::sequence ( void ) const {
  return string(_szInputSequence);
}

ListD PeptideFragment::a_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double dAion = 0.0,
      dBion = 0.0,
      dCion = 0.0,
      dXion = 0.0,
      dYion = 0.0,
      dZion = 0.0,
      dZdotion = 0.0,
      dNterm = pdMassAA['h'],
      // dCterm = pdMassAA['o'] + pdMassAA['h'],
      dProton = PROTON_MASS,
      dCO = pdMassAA['c'] + pdMassAA['o'],
      dH2 = pdMassAA['h'] + pdMassAA['h'],
      // dNH2 = pdMassAA['n'] + pdMassAA['h'] + pdMassAA['h'],
      dNH3 = pdMassAA['n'] + pdMassAA['h'] + pdMassAA['h'] + pdMassAA['h'];
   
    dBion = dNterm - pdMassAA['h'] + dProton;
    dYion = dPepMass;
   
    cerr << "OK" << endl;

    for ( i = 0; i < _iLenPeptide; i++ ) {
      cerr << "doing: " << i << "/" << _iLenPeptide << ": " << _szInputSequence[i] << endl;
      if ( i < _iLenPeptide - 1 ) {
        dBion += pdMassAA[(int)(_szInputSequence[i])];
        dBion += pdPositionMod[i];
        dAion = dBion - dCO;
        dCion = dBion + dNH3;
      }
      if (i > 0) {
        dYion -= pdMassAA[(int)(_szInputSequence[i-1])];
        dYion -= pdPositionMod[i-1];
        dXion = dYion + dCO - dH2;
        dZion = dYion - dNH3;
        dZdotion = dZion + pdMassAA['h'];
      }
      cerr << "OK" << endl;
   
      // stash A-ions
      if (i < _iLenPeptide - 1)
        a.push_back((dAion + (_iCharge-1)*dProton)/_iCharge);
      else
        a.push_back(0);  
   
      // stash B-ions
      if (i < _iLenPeptide - 1)
        b.push_back((dBion + (_iCharge-1)*dProton)/_iCharge);
      else
        b.push_back(0);  
   
      // stash C-ions
      if (i < _iLenPeptide - 1)
        c.push_back((dCion + (_iCharge-1)*dProton)/_iCharge);
      else
        c.push_back(0);
   
      // stash X-ions
      if ( i > 0 )
        x.push_back((dXion + (_iCharge-1)*dProton)/_iCharge);
      else
        x.push_back(0);

      // stash Y-ions
      if ( i > 0 )
        y.push_back((dYion + (_iCharge-1)*dProton)/_iCharge);
      else
        y.push_back(0);
   
      // stash Z-ions
      if (i > 0)
        z.push_back((dZion + (_iCharge-1)*dProton)/_iCharge);
      else
        z.push_back(0);

      // stash Z-dot-ions
      if ( i > 0 )
        zdot.push_back((dZdotion + (_iCharge-1)*dProton)/_iCharge);
      else
        zdot.push_back(0);
    }
    cerr << "OK again" << endl;
  }

  return a;
}


ListD PeptideFragment::b_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double dBion = 0.0;
    double dNterm = pdMassAA['h'];
    double dProton = PROTON_MASS;

    dBion = dNterm - pdMassAA['h'] + dProton;
   
    for ( i = 0; i < _iLenPeptide; i++ ) {
      if ( i < _iLenPeptide - 1 ) {
        dBion += pdMassAA[(int)(_szInputSequence[i])];
        dBion += pdPositionMod[i];
        b.push_back((dBion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        b.push_back(0);  
    }
  }

  return b;
}


ListD PeptideFragment::y_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double dYion = 0.0;
    double dProton = PROTON_MASS;
   
    dYion = dPepMass;
    for ( i = 0; i < _iLenPeptide; i++ ) {
      if (i > 0) {
        dYion -= pdMassAA[(int)(_szInputSequence[i-1])];
        dYion -= pdPositionMod[i-1];
        y.push_back((dYion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        y.push_back(0);
    }
  }

  return y;
}


double PeptideFragment::pI ( void ) {
  return compute_pI(_szInputSequence, _iLenPeptide, 0);
}


ListD PeptideFragment::peptide_mass ( void ) {
  double dProton = PROTON_MASS;
   
  mass.push_back(dPepMass - dProton);
  for ( int i = 1; i <= MAX_CHARGE; i++ )
    mass.push_back((dPepMass + (i-1)*dProton)/i);

  return mass;
}

string PeptideFragment::composition ( void ) {
  int iC, iH, iN, iO,iS;
  std::stringstream cat;

  iC=0;
  iH=0;
  iN=0;
  iO=0;
  iS=0;

  for ( unsigned int i = 0; i < _iLenPeptide; i++ ) {
    iC += piCompC[(unsigned int)_szInputSequence[i]];
    iH += piCompH[(unsigned int)_szInputSequence[i]];
    iN += piCompN[(unsigned int)_szInputSequence[i]];
    iO += piCompO[(unsigned int)_szInputSequence[i]];
    iS += piCompS[(unsigned int)_szInputSequence[i]];
  }

  cat << "C(" << iC << ") H(" << iH << ") N(" << iN << ") O(" << iO << ") S(" << iS << ")";
  return cat.str();
}

void PeptideFragment::_initialize_mass( void ) {
  double H, O, C, N, S;

  for (int i = 0; i < 128; i++)
    pdMassAA[i] = 999999.9;

  if (_iMassType) {
    H = pdMassAA['h'] =  1.0078250352; /* hydrogen */
    O = pdMassAA['o'] = 15.99491463;   /* oxygen */
    C = pdMassAA['c'] = 12.0000000;    /* carbon */
    N = pdMassAA['n'] = 14.003074;     /* nitrogen */
        pdMassAA['p'] = 30.973762;     /* phosphorus */
    S = pdMassAA['s'] = 31.9720707;    /* sulphur */
  }
  else { /* average masses */
    H = pdMassAA['h'] =  1.00794;  /* hydrogen */
    O = pdMassAA['o'] = 15.9994;   /* oxygen */
    C = pdMassAA['c'] = 12.0107;   /* carbon */
    N = pdMassAA['n'] = 14.0067;   /* nitrogen */
        pdMassAA['p'] = 30.973761; /* phosporus */
    S = pdMassAA['s'] = 32.065;    /* sulphur */
  }

  pdMassAA['G'] = C*2  + H*3  + N   + O;
  pdMassAA['A'] = C*3  + H*5  + N   + O;
  pdMassAA['S'] = C*3  + H*5  + N   + O*2;
  pdMassAA['P'] = C*5  + H*7  + N   + O;
  pdMassAA['V'] = C*5  + H*9  + N   + O;
  pdMassAA['T'] = C*4  + H*7  + N   + O*2;
  pdMassAA['C'] = C*3  + H*5  + N   + O   + S;
  pdMassAA['L'] = C*6  + H*11 + N   + O;
  pdMassAA['I'] = C*6  + H*11 + N   + O;
  pdMassAA['N'] = C*4  + H*6  + N*2 + O*2;
  pdMassAA['D'] = C*4  + H*5  + N   + O*3;
  pdMassAA['Q'] = C*5  + H*8  + N*2 + O*2;
  pdMassAA['K'] = C*6  + H*12 + N*2 + O;
  pdMassAA['E'] = C*5  + H*7  + N   + O*3;
  pdMassAA['M'] = C*5  + H*9  + N   + O   + S;
  pdMassAA['H'] = C*6  + H*7  + N*3 + O;
  pdMassAA['F'] = C*9  + H*9  + N   + O;
  pdMassAA['R'] = C*6  + H*12 + N*4 + O;
  pdMassAA['Y'] = C*9  + H*9  + N   + O*2;
  pdMassAA['W'] = C*11 + H*10 + N*2 + O;

  pdMassAA['O'] = C*5  + H*12 + N*2 + O*2;
  pdMassAA['X'] = pdMassAA['L'];  /* treat X as L or I for no good reason */
  pdMassAA['B'] = (pdMassAA['N'] + pdMassAA['D']) / 2.0;  /* treat B as average of N and D */
  pdMassAA['Z'] = (pdMassAA['Q'] + pdMassAA['E']) / 2.0;  /* treat Z as average of Q and E */
}


void PeptideFragment::_init_comp ( void ) {
  for (int i=0; i < 128; i++) {
    piCompC[i]=0;
    piCompH[i]=0;
    piCompN[i]=0;
    piCompO[i]=0;
    piCompS[i]=0;
  }

  piCompC['G'] = 2  ;
  piCompC['A'] = 3  ;
  piCompC['S'] = 3  ;
  piCompC['P'] = 5  ;
  piCompC['V'] = 5  ;
  piCompC['T'] = 4  ;
  piCompC['C'] = 3  ;
  piCompC['L'] = 6  ;
  piCompC['I'] = 6  ;
  piCompC['N'] = 4  ;
  piCompC['D'] = 4  ;
  piCompC['Q'] = 5  ;
  piCompC['K'] = 6  ;
  piCompC['E'] = 5  ;
  piCompC['M'] = 5  ;
  piCompC['H'] = 6  ;
  piCompC['F'] = 9  ;
  piCompC['R'] = 6  ;
  piCompC['Y'] = 9  ;
  piCompC['W'] = 11 ;
  piCompC['O'] = 5  ;

  piCompH['G'] = 3  ;
  piCompH['A'] = 5  ;
  piCompH['S'] = 5  ;
  piCompH['P'] = 7  ;
  piCompH['V'] = 9  ;
  piCompH['T'] = 7  ;
  piCompH['C'] = 5  ;
  piCompH['L'] = 11 ;
  piCompH['I'] = 11 ;
  piCompH['N'] = 6  ;
  piCompH['D'] = 5  ;
  piCompH['Q'] = 8  ;
  piCompH['K'] = 12 ;
  piCompH['E'] = 7  ;
  piCompH['M'] = 9  ;
  piCompH['H'] = 7  ;
  piCompH['F'] = 9  ;
  piCompH['R'] = 12 ;
  piCompH['Y'] = 9  ;
  piCompH['W'] = 10 ;
  piCompH['O'] = 12 ;

  piCompN['G'] = 1 ;
  piCompN['A'] = 1 ;
  piCompN['S'] = 1 ;
  piCompN['P'] = 1 ;
  piCompN['V'] = 1 ;
  piCompN['T'] = 1 ;
  piCompN['C'] = 1 ;
  piCompN['L'] = 1 ;
  piCompN['I'] = 1 ;
  piCompN['N'] = 2 ;
  piCompN['D'] = 1 ;
  piCompN['Q'] = 2 ;
  piCompN['K'] = 2 ;
  piCompN['E'] = 1 ;
  piCompN['M'] = 1 ;
  piCompN['H'] = 3 ;
  piCompN['F'] = 1 ;
  piCompN['R'] = 4 ;
  piCompN['Y'] = 1 ;
  piCompN['W'] = 2 ;
  piCompN['O'] = 2 ;

  piCompO['G'] = 1 ;
  piCompO['A'] = 1 ;
  piCompO['S'] = 2 ;
  piCompO['P'] = 1 ;
  piCompO['V'] = 1 ;
  piCompO['T'] = 2 ;
  piCompO['C'] = 1 ;
  piCompO['L'] = 1 ;
  piCompO['I'] = 1 ;
  piCompO['N'] = 2 ;
  piCompO['D'] = 3 ;
  piCompO['Q'] = 2 ;
  piCompO['K'] = 1 ;
  piCompO['E'] = 3 ;
  piCompO['M'] = 1 ;
  piCompO['H'] = 1 ;
  piCompO['F'] = 1 ;
  piCompO['R'] = 1 ;
  piCompO['Y'] = 2 ;
  piCompO['W'] = 1 ;
  piCompO['O'] = 2 ;

  piCompS['C'] = 1;
  piCompS['M'] = 1;
}
