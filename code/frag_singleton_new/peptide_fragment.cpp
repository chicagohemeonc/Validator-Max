#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include "peptide_fragment.hpp"
using namespace std;


PeptideFragment::PeptideFragment ( int massType ) {
  _iMassType = massType;
  _initialize_mass();
  _copy_mass();
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
  _copy_mass();

  /* remove non-alpha chars */
  j = 0;
  for (i = 0; i < strlen(peptide); i++)
    if (isalpha(peptide[i]))
      _szInputSequence[j++] = peptide[i];
  _szInputSequence[j] = '\0';
  _iLenPeptide = j;

  _iCharge = charge;
  _modification = string(modification);

  for (std::list<int>::iterator i = _modifiedPos.begin(); i != _modifiedPos.end(); ++i)
    pdPositionMod[*i] = 0;
  _modifiedPos.clear();

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

        sscanf(tok, "%lf %s", &dMass, szResidue); /* this is what remains after replacing '@' with ' ' */

        /*
         * first check if residue or number is entered
         */
        if (strspn(szResidue, "0123456789") == strlen(szResidue)) {
          int iPos=0;

          sscanf(szResidue, "%d", &iPos);
          if (iPos >= 1) {
            pdPositionMod[iPos-1] = dMass;
	    _modifiedPos.push_back(iPos-1);
	  }
        }
        else {
          pdMassAA[toupper(szResidue[0])] += dMass;
        }

        tok = strtok(NULL, " ");
      }
    }
  }

  dNterm = pdMassAA['h'] + pdMassAA['['];
  dCterm = pdMassAA['o'] + pdMassAA['h'] + pdMassAA[']'];

  dPepMass = dNterm + dCterm + dProton;
  for ( i = 0; i < _iLenPeptide; i++ ) {
    dPepMass += pdMassAA[(int)(_szInputSequence[i])] + pdPositionMod[i];
  }
}


std::string PeptideFragment::sequence ( void ) const {
  return string(_szInputSequence);
}

ListD PeptideFragment::a_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double
      dProton = PROTON_MASS,
      dAion = 0.0,
      dBion = dNterm - pdMassAA['h'] + dProton,
      dCO = pdMassAA['c'] + pdMassAA['o'];

    for ( i = 0; i < _iLenPeptide; i++ ) {
      if ( i < _iLenPeptide - 1 ) {
        dBion += pdMassAA[(int)(_szInputSequence[i])] + pdPositionMod[i];
        dAion = dBion - dCO;
        a.push_back((dAion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        a.push_back(0);
    }
  }

  return a;
}


ListD PeptideFragment::b_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double dBion = 0.0;
    double dProton = PROTON_MASS;

    dBion = dNterm - pdMassAA['h'] + dProton;

    for ( i = 0; i < _iLenPeptide; i++ ) {
      if ( i < _iLenPeptide - 1 ) {
        dBion += pdMassAA[(int)(_szInputSequence[i])] + pdPositionMod[i];
        b.push_back((dBion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        b.push_back(0);
    }
  }

  return b;
}


ListD PeptideFragment::c_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double dCion = 0.0;
    double dProton = PROTON_MASS;
    double dBion = dNterm - pdMassAA['h'] + dProton;
    double dNH3 = pdMassAA['n'] + pdMassAA['h'] + pdMassAA['h'] + pdMassAA['h'];

    for ( i = 0; i < _iLenPeptide; i++ ) {
      if ( i < _iLenPeptide - 1 ) {
        dBion += pdMassAA[(int)(_szInputSequence[i])] + pdPositionMod[i];
        dCion = dBion + dNH3;
        c.push_back((dCion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        c.push_back(0);
    }
  }

  return c;
}


ListD PeptideFragment::x_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double
      dProton = PROTON_MASS,
      dXion = 0.0,
      dYion = dPepMass,
      dCO = pdMassAA['c'] + pdMassAA['o'],
      dH2 = pdMassAA['h'] + pdMassAA['h'];

    for ( i = 0; i < _iLenPeptide; i++ ) {
      if ( i > 0 ) {
        dYion -= pdMassAA[(int)(_szInputSequence[i-1])] + pdPositionMod[i-1];
        dXion = dYion + dCO - dH2;
        x.push_back((dXion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        x.push_back(0);

      if ( i == 0 ) {
        dYion -= pdMassAA['['];
      }
    }
  }

  return x;
}


ListD PeptideFragment::y_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double dYion = dPepMass;
    double dProton = PROTON_MASS;

    for ( i = 0; i < _iLenPeptide; i++ ) {
      if (i > 0) {
        dYion -= pdMassAA[(int)(_szInputSequence[i-1])] + pdPositionMod[i-1];
        y.push_back((dYion + (_iCharge-1)*dProton)/_iCharge);
      }
      else {
        y.push_back(0);
      }

      if ( i == 0 ) {
        dYion -= pdMassAA['['];
      }
    }
  }

  return y;
}


ListD PeptideFragment::z_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double
      dProton = PROTON_MASS,
      dYion = dPepMass,
      dZion = 0.0,
      dNH3 = pdMassAA['n'] + pdMassAA['h'] + pdMassAA['h'] + pdMassAA['h'];

    for ( i = 0; i < _iLenPeptide; i++ ) {
      if ( i > 0 ) {
        dYion -= pdMassAA[(int)(_szInputSequence[i-1])] + pdPositionMod[i-1];
        dZion = dYion - dNH3;
        z.push_back((dZion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        z.push_back(0);

      if ( i == 0 )
	dYion -= pdMassAA['['];

    }
  }

  return z;
}


ListD PeptideFragment::zdot_ions ( void ) {

  if (_iLenPeptide > 0) {
    unsigned int i;
    double
      dProton = PROTON_MASS,
      dYion = dPepMass,
      dZion = 0.0,
      dZdotion = 0.0,
      dNH3 = pdMassAA['n'] + pdMassAA['h'] + pdMassAA['h'] + pdMassAA['h'];

    for ( i = 0; i < _iLenPeptide; i++ ) {
      if ( i > 0 ) {
        dYion -= pdMassAA[(int)(_szInputSequence[i-1])] + pdPositionMod[i-1];
        dZion = dYion - dNH3;
        dZdotion = dZion + pdMassAA['h'];
        zdot.push_back((dZdotion + (_iCharge-1)*dProton)/_iCharge);
      }
      else
        zdot.push_back(0);

      if ( i == 0 )
	dYion -= pdMassAA['['];

    }
  }

  return zdot;
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
    iC += piCompC[(int)_szInputSequence[i]];
    iH += piCompH[(int)_szInputSequence[i]];
    iN += piCompN[(int)_szInputSequence[i]];
    iO += piCompO[(int)_szInputSequence[i]];
    iS += piCompS[(int)_szInputSequence[i]];
  }

  cat << "C(" << iC << ") H(" << iH << ") N(" << iN << ") O(" << iO << ") S(" << iS << ")";
  return cat.str();
}

void PeptideFragment::_initialize_mass( void ) {
  double H=0.0, O=0.0, C=0.0, N=0.0, /* P=0.0, */ S=0.0;

  for (int i = 0; i < 128; i++)
    pdMassAA_template[i] = 999999.9;

  if (_iMassType) {
    H = pdMassAA_template['h'] =  1.0078250352; /* hydrogen */
    O = pdMassAA_template['o'] = 15.99491463;   /* oxygen */
    C = pdMassAA_template['c'] = 12.0000000;    /* carbon */
    N = pdMassAA_template['n'] = 14.003074;     /* nitrogen */
    // P = pdMassAA_template['p'] = 30.973762;     /* phosphorus */
    S = pdMassAA_template['s'] = 31.9720707;    /* sulphur */
  }
  else { /* average masses */
    H = pdMassAA_template['h'] =  1.00794;  /* hydrogen */
    O = pdMassAA_template['o'] = 15.9994;   /* oxygen */
    C = pdMassAA_template['c'] = 12.0107;   /* carbon */
    N = pdMassAA_template['n'] = 14.0067;   /* nitrogen */
    // P = pdMassAA_template['p'] = 30.973761; /* phosporus */
    S = pdMassAA_template['s'] = 32.065;    /* sulphur */
  }

  pdMassAA_template['G'] = C*2  + H*3  + N   + O;
  pdMassAA_template['A'] = C*3  + H*5  + N   + O;
  pdMassAA_template['S'] = C*3  + H*5  + N   + O*2;
  pdMassAA_template['P'] = C*5  + H*7  + N   + O;
  pdMassAA_template['V'] = C*5  + H*9  + N   + O;
  pdMassAA_template['T'] = C*4  + H*7  + N   + O*2;
  pdMassAA_template['C'] = C*3  + H*5  + N   + O   + S;
  pdMassAA_template['L'] = C*6  + H*11 + N   + O;
  pdMassAA_template['I'] = C*6  + H*11 + N   + O;
  pdMassAA_template['N'] = C*4  + H*6  + N*2 + O*2;
  pdMassAA_template['D'] = C*4  + H*5  + N   + O*3;
  pdMassAA_template['Q'] = C*5  + H*8  + N*2 + O*2;
  pdMassAA_template['K'] = C*6  + H*12 + N*2 + O;
  pdMassAA_template['E'] = C*5  + H*7  + N   + O*3;
  pdMassAA_template['M'] = C*5  + H*9  + N   + O   + S;
  pdMassAA_template['H'] = C*6  + H*7  + N*3 + O;
  pdMassAA_template['F'] = C*9  + H*9  + N   + O;
  pdMassAA_template['R'] = C*6  + H*12 + N*4 + O;
  pdMassAA_template['Y'] = C*9  + H*9  + N   + O*2;
  pdMassAA_template['W'] = C*11 + H*10 + N*2 + O;

  pdMassAA_template['O'] = C*5  + H*12 + N*2 + O*2;
  pdMassAA_template['X'] = pdMassAA_template['L'];  /* treat X as L or I for no good reason */
  pdMassAA_template['B'] = (pdMassAA_template['N'] + pdMassAA_template['D']) / 2.0;  /* treat B as average of N and D */
  pdMassAA_template['Z'] = (pdMassAA_template['Q'] + pdMassAA_template['E']) / 2.0;  /* treat Z as average of Q and E */

  pdMassAA_template['['] = 0.0;  /* initialize n-term mass */
  pdMassAA_template[']'] = 0.0;  /* initialize c-term mass */
}

void PeptideFragment::_copy_mass( void ) {
  memcpy ( pdMassAA, pdMassAA_template, 128*sizeof(double) );
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
