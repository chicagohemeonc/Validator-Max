#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include "peptide_fragment.hpp"
using namespace std;

const char* peptide[] = {
  "KMAEILVATVAFLLPSAEYSSVETDK",
  "KMAEILVATVAFLLPSAEYSSVETDKK",
  "MAEILVATVAFLLPSAEYSSVETDK",
  "MAEILVATVAFLLPSAEYSSVETDKK",
  "MAEILVATVAFLLPSAEYSSVETDKKFIVSLLLCLLDWCMALPVSVLLHPVSTAVLEEQHSAR",
  "KFIVSLLLCLLDWCMALPVSVLLHPVSTAVLEEQHSAR",
  "KFIVSLLLCLLDWCMALPVSVLLHPVSTAVLEEQHSARAPLLDYIYR",
  "FIVSLLLCLLDWCMALPVSVLLHPVSTAVLEEQHSAR",
  "FIVSLLLCLLDWCMALPVSVLLHPVSTAVLEEQHSARAPLLDYIYR",
  "FIVSLLLCLLDWCMALPVSVLLHPVSTAVLEEQHSARAPLLDYIYRVLHCCVCGSSTYTQQSHYILTLADLSSTDYDPFLPLANVK",
  "APLLDYIYR",
  "APLLDYIYRVLHCCVCGSSTYTQQSHYILTLADLSSTDYDPFLPLANVK",
  "APLLDYIYRVLHCCVCGSSTYTQQSHYILTLADLSSTDYDPFLPLANVKSSEPVQYHSSAELGNLLTVEEEK",
  "VLHCCVCGSSTYTQQSHYILTLADLSSTDYDPFLPLANVK",
  "VLHCCVCGSSTYTQQSHYILTLADLSSTDYDPFLPLANVKSSEPVQYHSSAELGNLLTVEEEK",
  "VLHCCVCGSSTYTQQSHYILTLADLSSTDYDPFLPLANVKSSEPVQYHSSAELGNLLTVEEEKK",
  "SSEPVQYHSSAELGNLLTVEEEK",
  "SSEPVQYHSSAELGNLLTVEEEKK",
  "SSEPVQYHSSAELGNLLTVEEEKKR",
  "RRSLELIPLTAR",
  "RSLELIPLTAR",
  "RSLELIPLTARMVMAHLVNHLGHYPLSGGPAILHSLVSENHDNAHVEGSELSFEVFR",
  "SLELIPLTAR",
  "SLELIPLTARMVMAHLVNHLGHYPLSGGPAILHSLVSENHDNAHVEGSELSFEVFR",
  "SLELIPLTARMVMAHLVNHLGHYPLSGGPAILHSLVSENHDNAHVEGSELSFEVFRSPNLQLFVFNDSTLISYLQTPTEGPVGGSPVGSLSDVR",
  "MVMAHLVNHLGHYPLSGGPAILHSLVSENHDNAHVEGSELSFEVFR",
  "MVMAHLVNHLGHYPLSGGPAILHSLVSENHDNAHVEGSELSFEVFRSPNLQLFVFNDSTLISYLQTPTEGPVGGSPVGSLSDVR",
  "MVMAHLVNHLGHYPLSGGPAILHSLVSENHDNAHVEGSELSFEVFRSPNLQLFVFNDSTLISYLQTPTEGPVGGSPVGSLSDVRVIVR",
  "SPNLQLFVFNDSTLISYLQTPTEGPVGGSPVGSLSDVR",
  "SPNLQLFVFNDSTLISYLQTPTEGPVGGSPVGSLSDVRVIVR",
  "SPNLQLFVFNDSTLISYLQTPTEGPVGGSPVGSLSDVRVIVRDISGK",
  "VIVRDISGK",
  "VIVRDISGKYSWDGK",
  "DISGK",
  "DISGKYSWDGK",
  "DISGKYSWDGKVLYGPLEGCLAPNGR",
  "YSWDGK",
  "YSWDGKVLYGPLEGCLAPNGR",
  "YSWDGKVLYGPLEGCLAPNGRNPSFLISSWHR",
  "VLYGPLEGCLAPNGR",
  "VLYGPLEGCLAPNGRNPSFLISSWHR",
  "VLYGPLEGCLAPNGRNPSFLISSWHRDTFGPQK",
  "NPSFLISSWHR",
  "NPSFLISSWHRDTFGPQK",
  "NPSFLISSWHRDTFGPQKDSSQVEEGDDVLDK",
  "DTFGPQK",
  "DTFGPQKDSSQVEEGDDVLDK",
  "DTFGPQKDSSQVEEGDDVLDKLLENIGHTSPECLLPSQLNLNEPSLTPCGMNYDQEK",
  "DSSQVEEGDDVLDK",
  "DSSQVEEGDDVLDKLLENIGHTSPECLLPSQLNLNEPSLTPCGMNYDQEK",
  "DSSQVEEGDDVLDKLLENIGHTSPECLLPSQLNLNEPSLTPCGMNYDQEKEIIEVILR",
  "LLENIGHTSPECLLPSQLNLNEPSLTPCGMNYDQEK",
  "LLENIGHTSPECLLPSQLNLNEPSLTPCGMNYDQEKEIIEVILR",
  "LLENIGHTSPECLLPSQLNLNEPSLTPCGMNYDQEKEIIEVILRQNAQEDEYIQSHNFDSAMK",
  "EIIEVILR",
  "EIIEVILRQNAQEDEYIQSHNFDSAMK",
  "EIIEVILRQNAQEDEYIQSHNFDSAMKVTSQGQPSPVEPR",
  "QNAQEDEYIQSHNFDSAMK",
  "QNAQEDEYIQSHNFDSAMKVTSQGQPSPVEPR",
  "QNAQEDEYIQSHNFDSAMKVTSQGQPSPVEPRGPFYFCR",
};

int main () {
  ListD mass;
  const size_t npeptides = sizeof(peptide) / sizeof(peptide[0]);

  PeptideFragment f(1);

  f.analyze((char *)"ACKRM", 3, (char *)"140@1 16.0@M 57.1@3");

  cout << "sequence: " << f.sequence() << endl;
  cout << "composition: " << f.composition() << endl;

  cout << "Peptide masses:";
  mass = f.peptide_mass();
  for (ListD::iterator i = mass.begin(); i != mass.end(); ++i) {
    cout << "  " << *i;
  }
  cout << endl;

  for (int n = 0; n < 1000; n++) {
  for (size_t i = 0; i < npeptides; ++i) {
//     cout << peptide[i] << endl;
    f.analyze(const_cast<char *>(peptide[i]), 0, const_cast<char *>(""));

//     cout << "sequence: " << f.sequence() << endl;
//     cout << "composition: " << f.composition() << endl;

    cout << "Peptide masses:";
    mass = f.peptide_mass();
    for (ListD::iterator i = mass.begin(); i != mass.end(); ++i) {
      cout << "  " << *i;
    }
    cout << endl;
  }
  }
  return 0;
}
