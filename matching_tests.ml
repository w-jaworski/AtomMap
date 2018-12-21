(*
 *  AtomMap: maps atoms in chemical reactions
 *  Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl> 
 *  Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw                

 *  Copyright (C) 2015-2017 Sara Szymkuc <saraszymkuc atSPAMfree gmail dot com>          
 *  Copyright (C) 2015-2017 Barbara Mikulak <basia dot mikulak atSPAMfree gmail dot com> 
 *  Copyright (C) 2015-2017 Institute of Organic Chemistry Polish Academy of Sciences     

 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *)

open Types


(*let _ = MatchingExec.map_atoms2 "200_OrganicSyntheses_Apr7_2015_" "results/200_OrganicSyntheses_Apr7_2015_mapped.tex" (Import.load_organic_syntheses2 ())
let _ = MatchingExec.map_atoms2 "przyklady_trudne_" "results/przyklady_trudne_mapped.tex" (Import.load_przyklady_trudne ())
let _ = MatchingExec.map_atoms2 "przyklady_trudne2_" "results/przyklady_trudne2_mapped.tex" (Import.load_przyklady_trudne2 ())
let _ = MatchingExec.map_atoms2 "Mapping_algorithm_train_set1_" "results/Mapping_algorithm_train_set1_mapped.tex" (Import.load_paths ())
let _ = MatchingExec.map_atoms2 "250_OrganicSyntheses_" "results/250_OrganicSyntheses_mapped.tex" (Import.load_organic_syntheses ())
let _ = MatchingExec.map_atoms2 "cukry_5_OHinProd_" "results/cukry_5_OHinProd_mapped.tex" (Import.load_reactions Import.w_5_ohinprod_filename)*)
(* let _ = MatchingExec.map_atoms2 "cukry_5_OHinSub_" "results/cukry_5_OHinSub_mapped.tex" (Import.load_reactions Import.w_5_ohinsub_filename) *)
(*let _ = MatchingExec.map_atoms2 "cukry_6_OHinProd_" "results/cukry_6_OHinProd_mapped.tex" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = MatchingExec.map_atoms2 "cukry_6_OHinSub_" "results/cukry_6_OHinSub_mapped.tex" (Import.load_reactions Import.w_6_ohinsub_filename)
let _ = MatchingExec.map_atoms2 "MetaCyc_" "results/MetaCyc_mapped.tex" (Import.load_metacyc ())
let _ = MatchingExec.map_atoms2 "RT_" "results/RT_mapped.tex" (Import.load_rt ())
let _ = MatchingExec.map_atoms2 "rxn_db_" "results/rxn_db_mapped.tex" (Import.load_rxn ())*)

(*let _ = MatchingExec.print_reactions "Mapping_algorithm_train_set1_" "results/Mapping_algorithm_train_set1.tex" (Import.load_paths ())
let _ = MatchingExec.print_reactions "przyklady_trudne_" "results/przyklady_trudne.tex" (Import.load_przyklady_trudne ())
let _ = MatchingExec.print_reactions "cukry_5_OHinProd_" "results/cukry_5_OHinProd_mapped.tex" (Import.load_reactions Import.w_5_ohinprod_filename)
let _ = MatchingExec.print_reactions "cukry_5_OHinSub_" "results/cukry_5_OHinSub_mapped.tex" (Import.load_reactions Import.w_5_ohinsub_filename)
let _ = MatchingExec.print_reactions "cukry_6_OHinProd_" "results/cukry_6_OHinProd_mapped.tex" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = MatchingExec.print_reactions "cukry_6_OHinSub_" "results/cukry_6_OHinSub_mapped.tex" (Import.load_reactions Import.w_6_ohinsub_filename)*)
(* let _ = print_reactions "results/OrganicSynthesesReactions.tex" (Import.load_organic_syntheses ())*)
(* let _ = print_reactions "200_Org_" "results/200_OrganicSyntheses_Apr7_2015.tex" (Import.load_organic_syntheses2 ()) *)
(* let _ = print_reactions "print_reactions_6_" "results/print_reactions_6.tex" (Import.load_przyklady_trudne ())  *)

let inv_stoi_H_reactions = [
  "ClC1CCCCC1=O.CO>>COC(=O)C1CCCC1.[Cl-]";
  "C[Si](C)(C)OC1=CCCC1.O=C1C=CC(=O)C=C1>>C[Si](C)(C)OC12CCCC1c1cc(O)ccc1O2";
  "[H][C@@]1(OC(=O)C(O)=C1O)[C@@H](O)CO.CC(=O)\\C=C/C=O>>[H][C@]12OC(=O)[C@](O)(c3ccc(C)o3)[C@@]1(O)OC[C@@H]2O.O";
(*   "CO[C@H]1C[C@H](S[Si](C)(C)C)[C@@H](CO1)OC(=O)N[Si](C)(C)C.CCc1cnc(O[Si](C)(C)C)nc1O[Si](C)(C)C>>CCC1=CN(C(CC2SC2CO)OC)C(=O)NC1=O.C[Si+](C)C.C[Si+](C)C.C[Si+](C)C.C[Si](C)(C)N.O=C=O"; *)
  "C1CC=CC1.ClC(=O)Cc1ccccc1>>O=C1[C@H]2CCC[C@H]2[C@H]1c1ccccc1.[Cl-]";
  "C=C1C=CC=CC1=C>>C1Cc2ccccc12";
  "NNc1ccccc1.CCC(C)=O>>CC1=C(C)C2=C(N1)C=CC=C2.N.O";
  ]

(*[{empty_record with reaction_smile = "COC(=O)C(C)(CCCC1(Cc2ccccc12)C#N)C=C>>[H][C@@]12CCc3ccccc3[C@@]1(CCCC2(C)C(=O)OC)C#N"}]*)
(* let _ = map_atoms2 "results/map_atoms_4.tex" (Xlist.map spurious2_reactions (fun s -> {empty_record with reaction_smile=s}))   *)
(*let _ = map_atoms2 "map_atoms_5_" "results/map_atoms_5.tex"
  [{empty_record with reaction_smile = "CCN(CC)CC.ClC(C1=CC=CC=C1I)=O.NC2=CC=CC=C2>>IC3=C(C(NC4=CC=CC=C4)=O)C=CC=C3.CC[NH+](CC)CC.[Cl-]"}]*)
(* let _ = map_atoms2 "map_atoms_7_" "results/map_atoms_7.tex" (Xlist.map inv_stoi_H_reactions (fun s -> {empty_record with reaction_smile=s})) *)
(*let _ = map_atoms2 "perycykliczne_" "results/perycykliczne.tex"
  [{empty_record with reaction_smile = "C\\C=C\\[C@H](O)C\\C=C\\c1ccccc1>>C[C@H](CC=O)[C@@H](C=C)c1ccccc1"}]*)

let trudne_wybrane_reactions = [
(*   "2","CC(=C)C=C.C[Si](C)(C)OC(=C)C(=O)c1ccccc1>>CC1=CCC(CC1)(O[Si](C)(C)C)C(=O)c1ccccc1";
  "29","CS(=O)c1ccc(Cl)cc1.FC(F)(F)C(=O)OC(=O)C(F)(F)F>>FC(F)(F)C(=O)OCSc1ccc(Cl)cc1.OC(=O)C(F)(F)F";*)
(*  "49","CC(CCC1=C(C)CCCC1(C)C)=CCO.OC=O>>C[C@]12CC[C@]3(C1)[C@](C)(CCCC3(C)C)C2OC=O.O";
   "50","CC(=C)[C@@H](O)[C@@H](O)C(C)=C.CC(C)=O>>CC(=C)[C@@H]1OC(C)(C)C[C@]1(C)C=O"; *)
(*   "51","C\\C=C(\\C)C1(C)OC(C)O[C@H]1C>>C[C@H]1O[C@@H](C)[C@](C)([C@H]1C)C(C)=O"; *)
  "26","CCCC[Sn](CCCC)(CCCC)CO[C@@H]1C[C@H](C=C1)N(CC=C)C(=O)OC(C)(C)C>>CC(C)(C)OC(=O)N(CC=C)[C@@H]1CC=C[C@H]1CO.CCCC[SnH](CCCC)CCCC";
  "34","OC1(C=CC(=O)C=C1)C#Cc1ccccc1.C=C1CC(=O)O1>>CC(=O)Cc1cc(ccc1O)C#Cc1ccccc1.O=C=O";
  "40","O=Cc1ccccc1.CC(=C)C[Si](C)(C)C.Cl>>CC(=C)CC(O)c1ccccc1.C[Si](C)(C)Cl";
  "41","O=Cc1ccccc1.C[Si](C)(C)CC=C>>OC(CC=C)c1ccccc1";
  "46","CC#CCCC1=CC=C(O1)[C@@H](O)CNC(C)=O>>[H][C@]1(CCC#CC)O[C@@]([H])(CNC(C)=O)C(=O)C=C1";
  ]
(* let _ = map_atoms2 "przyklady_trudne_wybrane_" "results/przyklady_trudne_wybrane.tex" (Xlist.map trudne_wybrane_reactions (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let trudne_pozostale_reactions = [
(*   "9","[H][C@@]1(OC(=O)C(O)=C1O)[C@@H](O)CO.CC(=O)\\C=C/C=O>>[H][C@]12OC(=O)[C@](O)(c3ccc(C)o3)[C@@]1(O)OC[C@@H]2O.O";  *)
(*   "17","O=C1C(=O)c2cccc3cccc1c23.CC(=O)CC(=O)CC(C)=O.C1C2C=CC1C=C2>>CC(=O)c1ccc(C(C)=O)c2-c3cccc4cccc(-c12)c34.[C-]#[O+].C1C=CC=C1.O.O"; *)
(*   "17b","CC(=O)C1=C2C3=CC=CC4=CC=CC(C2=C(C(C)=O)C1=O)=C34.C1C2C=CC1C=C2>>CC(=O)C1=CC=C(C(C)=O)C2=C1C1=CC=CC3=CC=CC2=C13.[C-]#[O+].C1C=CC=C1" *)
  ]
(* let _ = map_atoms2 "przyklady_trudne_pozostale_" "results/przyklady_trudne_pozostale.tex" (Xlist.map trudne_pozostale_reactions (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let trudne_nowe_reactions = [
(*   "7","CC1(O)CCCC1.CCCCN=[N+]=[N-]>>CCCCN1CCCCC1C.N#N"; *)
(*   "8","CC1CCCC(=O)C1C.CC(=O)C=C>>CC1CCCC2=CC(=O)CCC12C.O"; *)
(*   "10","NCC1CC1.ON=O>>C1CC=C1.N#N.O.O"; *)
(*   "38","OC(C(F)(Sc1ccccc1)F)=O.O=Cc2ccccc2.[C-]#[N+]c3ccccc3.Nc4ccccc4>>FC(Sc5ccccc5)(C(N(c6ccccc6)C(c7ccccc7)C(Nc8ccccc8)=O)=O)F.O"; *)
  "12","CO[C@H]1C[C@H](S[Si](C)(C)C)[C@@H](CO1)OC(=O)N[Si](C)(C)C.CCc1cnc(O[Si](C)(C)C)nc1O[Si](C)(C)C>>CCC1=CN(C(CC2SC2CO)OC)C(=O)NC1=O.C[Si+](C)C.C[Si+](C)C.C[Si+](C)C.C[Si](C)(C)N.O=C=O";
  ]
(* let _ = MatchingExec.map_atoms2 "przyklady_trudne_" "results/przyklady_trudne_nowe.tex" (Xlist.map trudne_nowe_reactions (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s}))  *)

let org_syn_niezmaczowane = [
  "5", "S[H].CC1=CC=C(N(=O)=O)C=C1>>NC2=CC=C(C=C2)C=O.O.[S]";
  "32", "ClP(Cl)(Cl)(Cl)Cl.NC(CC#N)=O>>N#CCC#N.ClP(Cl)(Cl)=O.[H]Cl.[H]Cl";
  "49", "[O+]#[C-].[O+]#[C-].O=C(C1=CC=CC(N(=O)=O)=C1/C=C/C2=CC=CO2)N(C)C>>CN(C(C3=CC=CC4=C3C=C(C5=CC=CO5)N4)=O)C.O=C=O.O=C=O";
  "50", "[C-]#[O+].[C-]#[O+].FC1=CC=C(C=C1)/C=C/C2=C(N(=O)=O)C=CC=C2C(N3CCCC3)=O>>FC4=CC=C(C5=CC6=C(C=CC=C6C(N7CCCC7)=O)N5)C=C4.O=C=O.O=C=O";
  "52", "[C-]#[O+].[C-]#[O+].O=C(C1=CC=CC(N(=O)=O)=C1/C=C/C2=CC3=C(C=C2)OCO3)N4CCOCC4>>O=C(C5=CC=CC6=C5C=C(C7=CC8=C(C=C7)OCO8)N6)N9CCOCC9.O=C=O.O=C=O";
  "59", "[C-]#[O+].[C-]#[O+].CS(=O)(N1CCN(C(C2=CC=CC(N(=O)=O)=C2/C=C/C3=CC(Cl)=CC(F)=C3)=O)CC1)=O>>CS(=O)(N4CCN(C(C5=CC=CC6=C5C=C(C7=CC(Cl)=CC(F)=C7)N6)=O)CC4)=O.O=C=O.O=C=O";
  "61", "[C-]#[O+].[C-]#[O+].CC(C)(OC(N1CCN(C(C2=CC=CC(N(=O)=O)=C2/C=C/C3=CC4=C(C=C3)OCO4)=O)CC1)=O)C>>CC(C)(OC(N5CCN(C(C6=CC=CC7=C6C=C(C8=CC9=C(C=C8)OCO9)N7)=O)CC5)=O)C.O=C=O.O=C=O";
  "80", "[C-]#[O+].[C-]#[O+].FC1=CC=C(C=C1)/C=C/C2=C(N(=O)=O)C=CC=C2C(N3CCC(N4C(NC5=C4C=CC=C5)=O)CC3)=O>>FC6=CC=C(C7=CC8=C(C=CC=C8C(N9CCC(N%10C(NC%11=C%10C=CC=C%11)=O)CC9)=O)N7)C=C6.O=C=O.O=C=O";
  "88", "O=Cc1ccccc1.FS(F)(c2ccccc2)F>>FC(c3ccccc3)F.FS(c4ccccc4)=O";
  "105", "ClP(Cl)(Cl)=O.O=C1C2=C(NC3=C1C=CC=C3)C=CC=C2>>ClC4=C5C=CC=CC5=NC6=C4C=CC=C6.OP(Cl)(Cl)=O";
  "125", "O=C1OC(C=C1)=O.NC2=CC=CC=C2>>OC(/C=C\\C(NC3=CC=CC=C3)=O)=O";
  "164", "BrC1=CC=CC=C1.O=C(CC2)OC2=O>>BrC3=CC=C(C(CCC(O)=O)=O)C=C3";
  "171", "OCCOCCN.CC(OC(C)=O)=O>>CC(O)=O.OCCOCCNC(C)=O";
  ]
(* let _ = MatchingExec.map_atoms2 "200_OrganicSyntheses_Apr7_2015_" "results/org_syn_niezmaczowane.tex" (Xlist.map org_syn_niezmaczowane (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let org_syn_estry = [
   "6", "OCCCCCO.CC(OC(C)=O)=O.CC(OC(C)=O)=O>>CC(OCCCCCOC(C)=O)=O.CC(O)=O.CC(O)=O";
  "91", "CCOC(C1CC(C(C(OCC)=O)CC1=O)=O)=O.O.O>>O=C(CC2)CCC2=O.CCO.CCO.O=C=O.O=C=O";
  "125", "O=C1OC(C=C1)=O.NC2=CC=CC=C2>>OC(/C=C\\C(NC3=CC=CC=C3)=O)=O";
  "156", "O=C1OC(C2CC=CCC21)=O.OCC.CCO>>O=C(OCC)C3C(C(OCC)=O)CC=CC3.O";
  "164", "BrC1=CC=CC=C1.O=C(CC2)OC2=O>>BrC3=CC=C(C(CCC(O)=O)=O)C=C3";
   "171", "OCCOCCN.CC(OC(C)=O)=O>>CC(O)=O.OCCOCCNC(C)=O";
  ]
(* let _ = MatchingExec.map_atoms2 "200_OrganicSyntheses_Apr7_2015_" "results/org_syn_estry.tex" (Xlist.map org_syn_estry (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s}))  *)

let trudne2_2 = [
  "2", "COC(=O)C1C=CC1C(=O)OC>>COC(=O)\\C=C\\C=C/C(=O)OC";
  ]
(* let _ = MatchingExec.map_atoms2 "przyklady_trudne2_" "results/przyklady_trudne2_mapped.tex" (Xlist.map trudne2_2 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let trudne2_1 = [
  "1", "CN1C=C(C)C(=O)NC1=O>>[H][C@]1(N(C)C(=O)NC(=O)[C@]1(C)O)C1=NC(=O)N(C)C=C1C";
  ]
(* let _ = MatchingExec.map_atoms2 "przyklady_trudne2_" "results/przyklady_trudne2_1_mapped.tex" (Xlist.map trudne2_1 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let org_syn_niezmaczowane2 = [
(*  "32", "ClP(Cl)(Cl)(Cl)Cl.NC(CC#N)=O>>N#CCC#N.ClP(Cl)(Cl)=O.[H]Cl.[H]Cl";*)
(*   "88", "O=Cc1ccccc1.FS(F)(c2ccccc2)F>>FC(c3ccccc3)F.FS(c4ccccc4)=O"; *)
(*   "105", "ClP(Cl)(Cl)=O.O=C1C2=C(NC3=C1C=CC=C3)C=CC=C2>>ClC4=C5C=CC=CC5=NC6=C4C=CC=C6.OP(Cl)(Cl)=O"; *)
(*    "19", "COC(C(C)C)=O.FC1=CC(Br)=CC=C1.C2(CCCCC2)[N-]C3CCCCC3>>COC(C(C)(C4=CC(F)=CC=C4)C)=O.C5(CCCCC5)NC6CCCCC6.[Br-]";  *)
(*   "94", "CC(OC(C)=O)=O.CC(OC(C)=O)=O.O=Cc1ccccc1.OC(CCC(c2ccccc2)=O)=O>>O=C3OC(c4ccccc4)=C/C3=C\\c5ccccc5.CC(O)=O.CC(O)=O.CC(O)=O.CC(O)=O"; *) (* timeout *)
"195", "O=C(CC1=CC=CC=C1)CC2=CC=CC=C2.OC(C)=O.OC(C)=O>>O=C3C(C4=CC=CC=C4)=C(C)OC(C)=C3C5=CC=CC=C5.O.O.O";

(*  "6", "OCCCCCO.CC(OC(C)=O)=O.CC(OC(C)=O)=O>>CC(OCCCCCOC(C)=O)=O.CC(O)=O.CC(O)=O"; (* prawdopodobnie błąd przy redukcji rozwiązań symetrycznych *)
  "36", "CS(=O)(Cl)=O.CCN(CC)CC.CC(C)([C@H](NC(C1=C(C=CC=C1)Br)=O)CO)C>>CC(C)([C@H]2COC(C3=C(C=CC=C3)Br)=N2)C.CS(=O)(O)=O.CC[NH+](CC)CC.[Cl-]";
  "57", "CC([O-])=O.CC1(C)O[I](Cl)C2=CC=CC=C12.C[Si](C(F)(F)F)(C)C>>CC3(C)O[I](C4=C3C=CC=C4)C(F)(F)F.C[Si](OC(C)=O)(C)C.[Cl-]";
  "91", "CCOC(C1CC(C(C(OCC)=O)CC1=O)=O)=O.O.O>>O=C(CC2)CCC2=O.CCO.CCO.O=C=O.O=C=O";
  "119", "ClS(Cl)=O.NCCC1=C(C=C(C=C1)Cl)CO.[OH-].[OH-]>>ClC2=CC=C3CCNCC3=C2.O=S=O.O.O.[Cl-].[Cl-]";
  "120", "NC(CO)C1=CC=CC=C1.ClS(Cl)=O>>[NH3+]C(CCl)C2=CC=CC=C2.O=S=O.[Cl-]";
  "125", "O=C1OC(C=C1)=O.NC2=CC=CC=C2>>OC(/C=C\\C(NC3=CC=CC=C3)=O)=O";
  "126", "O=C(OC(C)=O)C.OC(/C=C\\C(NC1=CC=CC=C1)=O)=O>>O=C2N(C3=CC=CC=C3)C(C=C2)=O.O=C(O)C.O=C(O)C";
  "141", "O=C(OC)C1OC1(C)CC(OC)OC>>O=C(OC)C2=C(C)C=CO2.OC.OC";
  "149", "OCC1=CC=CO1.O=C(C2=CC=CC(Cl)=C2)OO>>O=C3COC(C=C3)O.OC(C4=CC=CC(Cl)=C4)=O";
  "151", "C1=CC=CC1.C2C=C2>>C34C=CC(C5C4C5)C3";
  "156", "O=C1OC(C2CC=CCC21)=O.OCC.CCO>>O=C(OCC)C3C(C(OCC)=O)CC=CC3.O";
  "160", "C=C/C=C/OC(C)=O.C1=CC=CCC=C1>>[H][C@@]23C=CC=C[C@@](CC=C[C@@H]3OC(C)=O)([H])C2";
  "164", "BrC1=CC=CC=C1.O=C(CC2)OC2=O>>BrC3=CC=C(C(CCC(O)=O)=O)C=C3";
  "169", "C[Si](C)(C)O/C=C1CCCCCCCCCCC1(CC(C)=C)O>>CC2=CC=C3CCCCCCCCCCC3=C2.O[Si](C)(C)C.O";
  "171", "OCCOCCN.CC(OC(C)=O)=O>>CC(O)=O.OCCOCCNC(C)=O";
  "193", "c1ccccc1CC#N.O=C(OC)OC>>CC(C#N)C1=CC=CC=C1.CO.O=C=O";
  "198", "CCCCCC/C=C\\CCCCCC(OCC1OC(C)(C)OC1)=O.O>>CCCCCC/C=C\\CCCCCC(OCC(CO)O)=O.CC(C)=O";  (* nazwa przegrupowania? *)
  *)
  ]
(* let _ = MatchingExec.map_atoms2 "200_OrganicSyntheses_Apr7_2015_" "results/org_syn_niezmaczowane2d.tex" (Xlist.map org_syn_niezmaczowane2 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let trudne2_poprawione = [
(* Stack overflow *)
(*   "20", "CCOC(=O)CCCCC(=O)OCC>>O=C1CCCC1"; *)
  "42", "[H][C@@]12CC(=O)O[C@@]3(C)[C@]1([H])[C@@]2(CC[C@@]3(C)O[Si](C)(C)C(C)(C)C)C1=CC=CCC1>>[H][C@]12C=CCCC1=C1CC[C@@]([H])(O[Si](C)(C)C(C)(C)C)[C@@]3(C)COC(=O)[C@]2([H])[C@@]13[H]";
(*  "48", "CCOC=C.CC#CC1(CC1)C(C)=O.O>>CC(=O)C1=CC[C@]2(C)CC[C@]12O.CCO";*)
(*   manage\_aromaticity *)
(*   "13", "Coc1ccc(cc1)C#Cc1ccc(OC)cc1.C=C>>COc1ccc(cc1)C(=C)C(=C)c1ccc(OC)cc1";  *)
(*  "14", "CC(=O)C=C.CCCCC=C>>CCCCC=CC(C)=O.C=C";
  "15", "[H]Oc1ccccc1.CCOC(=O)CC(=O)CC>>CCC1=CC(=O)Oc2ccccc12.CCO";
  "30", "CC[C@@H]1CCCN(C1C=C)C(=O)CC>>CC[C@@H]1CCCNC(=O)[C@H](C)C\\C=C1";
  "41", "C[N+]1([O-])CCCC1C1=CC=CN=C1.CC(=O)OC(C)=O>>CC(=O)N1CCCC1C1=CC=CN=C1.C=O.CC(O)=O";*)
(*  "43", "[H][C@]12CC[C@@]3(C)[C@@H](CC[C@@]3(O)C1COC(C2)OC)OS(=O)(=O)C1=CC=C(C)C=C1>>[H][C@@]12CC(OC)OCC1C(=O)CC\\C=C(C)\\CC2.CC1=CC=C(C=C1)S(O)(=O)=O";
  "47", "CC(=O)C1(CC1)C#Cc1ccccc1>>C[C@@]12CC[C@]1(O)C(=CC2)C(=O)c1ccccc1";
  "52", "C1N=C1c1ccccc1.O=C1C(C(C(C1c1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1>>c1ccc(cc1)C1C(=NC(c2ccccc2)=C(c2ccccc2)C(c2ccccc2)=C1c1ccccc1)c1ccccc1";
  "57", "O=C1C(=C(C(=C1c1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1.N#Cc1ccccc1>>c1ccc(cc1)-c1nc(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c1-c1ccccc1.[C-]#[O+]";
  "58", "O=C1C(=C(C(=C1c1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1.N#Cc1ccccc1>>c1ccc(cc1)-c1nc(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c1-c1ccccc1";*)
  ]
(* let _ = MatchingExec.map_atoms2 "przyklady_trudne2_" "results/przyklady_trudne2_mapped.tex" (Xlist.map trudne2_poprawione (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s}))  *)


(*
./map_reactions 60 "przyklady_trudne_" "../src_data/przyklady_trudne.dat" "results/przyklady_trudne_mapped.tex"
./map_reactions 60 "200_OrganicSyntheses_Apr7_2015_" ../src_data/200_OrganicSyntheses_Apr7_2015.dat "results/200_OrganicSyntheses_Apr7_2015_mapped.tex"
./map_reactions 60 "przyklady_trudne2_" "../src_data/reakcje_do_zmaczowania_trudne_2.dat" "results/przyklady_trudne2_mapped.tex"
./map_reactions 60 "Mapping_algorithm_train_set1_" "../src_data/Mapping_algorithm_train_set1.dat" "results/Mapping_algorithm_train_set1_mapped.tex"

./map_reactions 60 "przyklady_trudne_" "../src_data/przyklady_trudne.dat" "results/przyklady_trudne_mapped.tex" >log1.txt &
./map_reactions 60 "200_OrganicSyntheses_Apr7_2015_" ../src_data/200_OrganicSyntheses_Apr7_2015.dat "results/200_OrganicSyntheses_Apr7_2015_mapped.tex" >log2.txt &
./map_reactions 60 "przyklady_trudne2_" "../src_data/reakcje_do_zmaczowania_trudne_2.dat" "results/przyklady_trudne2_mapped.tex" >log3.txt &
./map_reactions 60 "Mapping_algorithm_train_set1_" "../src_data/Mapping_algorithm_train_set1.dat" "results/Mapping_algorithm_train_set1_mapped.tex" >log4.txt &

./map_reactions 60 "W_5_OHinProd_" "../src_data/W_5_OHinProd.txt" "results/W_5_OHinProd.tex"
./map_reactions 60 "W_5_OHinSub_" "../src_data/W_5_OHinSub.txt" "results/W_5_OHinSub.tex"
./map_reactions 60 "W_6_OHinProd_" "../src_data/W_6_OHinProd.txt" "results/W_6_OHinProd.tex"
./map_reactions 60 "W_6_OHinSub_" "../src_data/W_6_OHinSub.txt" "results/W_6_OHinSub.tex"
*)

let print_errors filename errors =
(*   if errors = [] then () else *)
  File.file_out filename (fun file ->
    Xlist.iter (List.rev errors) (fun (s,re) ->
      Printf.fprintf file "%s %s\n   %s\n" re.rxn_id re.reaction_smile s))

let print_results filename results =
  let results = Xlist.fold results [] (fun results (t,re) -> if t < 0.01 then results else (t,re) :: results) in
  let results = Xlist.sort results (fun (t1,_) (t2,_) -> compare t1 t2) in
  File.file_out filename (fun file ->
    Xlist.iter results (fun (t,re) ->
      Printf.fprintf file "%f %s %s\n" t re.rxn_id re.reaction_smile))

let print_results2 filename results =
(*   let results = Xlist.fold results [] (fun results (t,s,re) -> if t < 0.01 then results else (t,re) :: results) in *)
  let results = Xlist.sort results (fun (t1,_,_) (t2,_,_) -> compare t1 t2) in
  File.file_out filename (fun file ->
    Xlist.iter results (fun (t,size,re) ->
      Printf.fprintf file "%f %d %s %s\n" t size re.rxn_id re.reaction_smile))

let print_results3 filename results =
(*   let results = Xlist.fold results [] (fun results (t,s,re) -> if t < 0.01 then results else (t,re) :: results) in *)
  let results = Xlist.sort results (fun (t1,_,_,_,_) (t2,_,_,_,_) -> compare t1 t2) in
  File.file_out filename (fun file ->
    Xlist.iter results (fun (t,alt_labels_size,no_broken_bonds,rgood_cands_size,re) ->
      Printf.fprintf file "%f alt=%d no_br=%d gcand=%d %s %s\n" t alt_labels_size no_broken_bonds rgood_cands_size re.rxn_id re.reaction_smile))

let print_results4 filename results =
(*   let results = Xlist.fold results [] (fun results (t,s,re) -> if t < 0.01 then results else (t,re) :: results) in *)
  let results = Xlist.sort results (fun (t1,_,_,_,_) (t2,_,_,_,_) -> compare t1 t2) in
  File.file_out filename (fun file ->
    Xlist.iter results (fun (t,no_sols,quality,no_broken_bonds,re) ->
      Printf.fprintf file "%f sol=%d qual=%d no_br=%d %s %s\n" t no_sols quality no_broken_bonds re.rxn_id re.reaction_smile))

let test_prepare_record_for_matching res_filename err_filename reactions =
  let results,errors = Xlist.fold reactions ([],[]) (fun (results,errors) re ->
    try
      let time1 = Sys.time () in
      let (*r,messages*)_ = Smiles.prepare_record_for_matching re [] in
      let time2 = Sys.time () in
      (time2 -. time1,re) :: results, errors
    with e -> results, (Printexc.to_string e,re) :: errors) in
  print_results res_filename results;
  print_errors err_filename errors;
  ()

(*let _ = test_prepare_record_for_matching "results/results1_200_OrganicSyntheses_Apr7_2015.txt" "results/errors_200_OrganicSyntheses_Apr7_2015.txt" (Import.load_organic_syntheses2 ())
let _ = test_prepare_record_for_matching "results/results1_przyklady_trudne.txt" "results/errors_przyklady_trudne.txt" (Import.load_przyklady_trudne ())
let _ = test_prepare_record_for_matching "results/results1_przyklady_trudne2.txt" "results/errors_przyklady_trudne2.txt" (Import.load_przyklady_trudne2 ())
let _ = test_prepare_record_for_matching "results/results1_Mapping_algorithm_train_set1.txt" "results/errors_Mapping_algorithm_train_set1.txt" (Import.load_paths ())
let _ = test_prepare_record_for_matching "results/results1_250_OrganicSyntheses.txt" "results/errors_250_OrganicSyntheses.txt" (Import.load_organic_syntheses ())
let _ = test_prepare_record_for_matching "results/results1_cukry_5_OHinProd.txt" "results/errors_cukry_5_OHinProd.txt" (Import.load_reactions Import.w_5_ohinprod_filename)
let _ = test_prepare_record_for_matching "results/results1_cukry_5_OHinSub.txt" "results/errors_cukry_5_OHinSub.txt" (Import.load_reactions Import.w_5_ohinsub_filename)
let _ = test_prepare_record_for_matching "results/results1_cukry_6_OHinProd.txt" "results/errors_cukry_6_OHinProd.txt" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = test_prepare_record_for_matching "results/results1_cukry_6_OHinSub.txt" "results/errors_cukry_6_OHinSub.txt" (Import.load_reactions Import.w_6_ohinsub_filename)*)
(* let _ = test_prepare_record_for_matching "results/results1_MetaCyc.txt" "results/errors_MetaCyc.txt" (Import.load_metacyc ()) *)
(* let _ = test_prepare_record_for_matching "results/results1_RT.txt" "results/errors_RT.txt" (Import.load_rt ()) *)
(*let _ = test_prepare_record_for_matching "results/results1_rxn_db.txt" "results/errors_rxn_db.txt" (Import.load_rxn ())*)

let test_multilevel_label_reaction3 res_filename err_filename reactions =
  let results,errors = Xlist.fold reactions ([],[]) (fun (results,errors) re ->
    try
      let time1 = Sys.time () in
      Types.time := time1;
      let r,messages = Smiles.prepare_record_for_matching re [] in
      let r = MatchingExec.set_unbreakable r in
      let result,history,alt_labels_list = CommonSubstructure.multilevel_label_reaction3 r [] [] in
      let time2 = Sys.time () in
      (time2 -. time1,Xlist.size alt_labels_list,re) :: results, errors
    with e -> results, (Printexc.to_string e,re) :: errors) in
  print_results2 res_filename results;
  print_errors err_filename errors;
  ()

(*let _ = test_multilevel_label_reaction3 "results/results2_200_OrganicSyntheses_Apr7_2015.txt" "results/errors_200_OrganicSyntheses_Apr7_2015.txt" (Import.load_organic_syntheses2 ())
let _ = test_multilevel_label_reaction3 "results/results2_przyklady_trudne.txt" "results/errors_przyklady_trudne.txt" (Import.load_przyklady_trudne ())
let _ = test_multilevel_label_reaction3 "results/results2_przyklady_trudne2.txt" "results/errors_przyklady_trudne2.txt" (Import.load_przyklady_trudne2 ())
let _ = test_multilevel_label_reaction3 "results/results2_Mapping_algorithm_train_set1.txt" "results/errors_Mapping_algorithm_train_set1.txt" (Import.load_paths ())
let _ = test_multilevel_label_reaction3 "results/results2_250_OrganicSyntheses.txt" "results/errors_250_OrganicSyntheses.txt" (Import.load_organic_syntheses ())
let _ = test_multilevel_label_reaction3 "results/results2_cukry_5_OHinProd.txt" "results/errors_cukry_5_OHinProd.txt" (Import.load_reactions Import.w_5_ohinprod_filename)
let _ = test_multilevel_label_reaction3 "results/results2_cukry_5_OHinSub.txt" "results/errors_cukry_5_OHinSub.txt" (Import.load_reactions Import.w_5_ohinsub_filename)
let _ = test_multilevel_label_reaction3 "results/results2_cukry_6_OHinProd.txt" "results/errors_cukry_6_OHinProd.txt" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = test_multilevel_label_reaction3 "results/results2_cukry_6_OHinSub.txt" "results/errors_cukry_6_OHinSub.txt" (Import.load_reactions Import.w_6_ohinsub_filename)
let _ = test_multilevel_label_reaction3 "results/results2_MetaCyc.txt" "results/errors_MetaCyc.txt" (Import.load_metacyc ())
let _ = test_multilevel_label_reaction3 "results/results2_RT.txt" "results/errors_RT.txt" (Import.load_rt ())
let _ = test_multilevel_label_reaction3 "results/results2_rxn_db.txt" "results/errors_rxn_db.txt" (Import.load_rxn ())*)

(* let _ = MatchingExec.print_labeled_reactions "Mapping_algorithm_train_set1_" "results/labeled_Mapping_algorithm_train_set1.tex" (Import.load_paths ()) *)
(* let _ = MatchingExec.print_labeled_reactions "cukry_5_OHinProd_" "results/labeled_cukry_5_OHinProd.tex" (Import.load_reactions Import.w_5_ohinprod_filename) *)
(* let _ = MatchingExec.print_labeled_reactions "cukry_5_OHinSub_" "results/labeled_cukry_5_OHinSub.tex" (Import.load_reactions Import.w_5_ohinsub_filename) *)
(*let _ = MatchingExec.print_labeled_reactions "cukry_6_OHinProd_" "results/labeled_cukry_6_OHinProd.tex" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = MatchingExec.print_labeled_reactions "cukry_6_OHinSub_" "results/labeled_cukry_6_OHinSub.tex" (Import.load_reactions Import.w_6_ohinsub_filename)*)

let test_create_good_candidates2 res_filename err_filename reactions =
  let results,errors = Xlist.fold reactions ([],[]) (fun (results,errors) re ->
    try
      Printf.printf "%s %s\n%!" re.rxn_id re.reaction_smile;
      let time1 = Sys.time () in
      Types.time := time1;
      let r,messages = Smiles.prepare_record_for_matching re [] in
      let r = MatchingExec.set_unbreakable r in
      let msg = ref messages in
      let result,history,alt_labels_list = CommonSubstructure.multilevel_label_reaction3 r [] [] in
      let alt_labels_size = Xlist.size alt_labels_list in
      let alt_labels_list = if alt_labels_list = [] then [r.empty_labels] else alt_labels_list in
      let arl = [r,Collection.of_list alt_labels_list] in
      let no_broken_bonds, _, rgood_cands = AtomMapping.create_good_candidates2 msg 6 0 arl in
      let time2 = Sys.time () in
      (time2 -. time1,alt_labels_size,no_broken_bonds,Xlist.size rgood_cands,re) :: results, errors
    with e -> results, (Printexc.to_string e,re) :: errors) in
  print_results3 res_filename results;
  print_errors err_filename errors;
  ()

(*let _ = test_create_good_candidates2 "results/results3_200_OrganicSyntheses_Apr7_2015.txt" "results/errors_200_OrganicSyntheses_Apr7_2015.txt" (Import.load_organic_syntheses2 ())
let _ = test_create_good_candidates2 "results/results3_przyklady_trudne.txt" "results/errors_przyklady_trudne.txt" (Import.load_przyklady_trudne ())
let _ = test_create_good_candidates2 "results/results3_przyklady_trudne2.txt" "results/errors_przyklady_trudne2.txt" (Import.load_przyklady_trudne2 ())*)
(* let _ = test_create_good_candidates2 "results/results3_Mapping_algorithm_train_set1.txt" "results/errors_Mapping_algorithm_train_set1.txt" (Import.load_paths ()) *)
(*let _ = test_create_good_candidates2 "results/results3_250_OrganicSyntheses.txt" "results/errors_250_OrganicSyntheses.txt" (Import.load_organic_syntheses ())*)
(*let _ = test_create_good_candidates2 "results/results3_cukry_5_OHinProd.txt" "results/errors_cukry_5_OHinProd.txt" (Import.load_reactions Import.w_5_ohinprod_filename)*)
(* let _ = test_create_good_candidates2 "results/results3_cukry_5_OHinSub.txt" "results/errors_cukry_5_OHinSub.txt" (Import.load_reactions Import.w_5_ohinsub_filename)  *)
(*let _ = test_create_good_candidates2 "results/results3_cukry_6_OHinProd.txt" "results/errors_cukry_6_OHinProd.txt" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = test_create_good_candidates2 "results/results3_cukry_6_OHinSub.txt" "results/errors_cukry_6_OHinSub.txt" (Import.load_reactions Import.w_6_ohinsub_filename)*)
(*let _ = test_create_good_candidates2 "results/results3_MetaCyc.txt" "results/errors_MetaCyc.txt" (Import.load_metacyc ())
let _ = test_create_good_candidates2 "results/results3_RT.txt" "results/errors_RT.txt" (Import.load_rt ())
let _ = test_create_good_candidates2 "results/results3_rxn_db.txt" "results/errors_rxn_db.txt" (Import.load_rxn ())*)

let pamieciozerne = [
(*   "3", "CCCC(=O)O[C@H]1[C@@H]([C@H]2COC(C)(C)O2)O[C@@H]2OC(C)(C)O[C@@H]21>>CCCC(=O)O[C@@H]1[C@@H](O)C(O)O[C@@H]1[C@H]1COC(C)(C)O1.CCCC(=O)O[C@H]1[C@@H]([C@H](O)CO)O[C@@H]2OC(C)(C)O[C@@H]21" *)
(*   "15", "COC1(CO)OC(CO)C(O)C1O.CC1OC(c2ccccc2)OC2C(O)C(CO[Si](C)(C)C(C)(C)C)OC21O>>COC12COC(c3ccccc3)OC1C(O)C(CO)O2" *)
(*   "1", "CC1(C)OC[C@H]([C@H]2OC(O)[C@H]3OC(C)(C)O[C@@H]23)O1.CO[C@@H]1O[C@H](COS(=O)(=O)C(F)(F)F)[C@H]2OC(C)(C)O[C@@H]12>>CO[C@@H]1O[C@H](CO[C@H]2O[C@H]([C@H]3COC(C)(C)O3)[C@@H]3OC(C)(C)O[C@H]23)[C@H]2OC(C)(C)O[C@@H]12" *)
  "3427", "CC(C)(C)c1cc(CCC(CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)NCCNC(CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)cc(C(C)(C)C)c1O>>CC(=O)N(CCN(C(C)=O)C(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)C(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1";
  "3428", "CC(CNC(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)NC(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1>>CC(=O)N(CC(C)N(C(C)=O)C(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)C(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1";
  "3429", "CC(C)(C)c1cc(CCC(CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)NCCCCCCNC(CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)cc(C(C)(C)C)c1O>>CC(=O)N(CCCCCCN(C(C)=O)C(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)C(CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1";
  ]

(* let _ = MatchingExec.print_mapped_reactions "Mapping_algorithm_train_set1_" "results/mapped_Mapping_algorithm_train_set1.tex" (Import.load_paths ()) *)
(*let _ = MatchingExec.print_mapped_reactions "cukry_5_OHinProd_" "results/mapped_cukry_5_OHinProd.tex" (Import.load_reactions Import.w_5_ohinprod_filename)
let _ = MatchingExec.print_mapped_reactions "cukry_5_OHinSub_" "results/mapped_cukry_5_OHinSub.tex" (Import.load_reactions Import.w_5_ohinsub_filename)*)
(*let _ = MatchingExec.print_mapped_reactions "cukry_6_OHinProd_" "results/mapped_cukry_6_OHinProd.tex" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = MatchingExec.print_mapped_reactions "cukry_6_OHinSub_" "results/mapped_cukry_6_OHinSub.tex" (Import.load_reactions Import.w_6_ohinsub_filename) *)
(* let _ = MatchingExec.print_mapped_reactions "alk_" "results/alk.tex" (Xlist.map pamieciozerne (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let test_broken_bonds n res_filename err_filename reactions =
  let results,errors,_ = Xlist.fold reactions ([],[],n) (fun (results,errors,n) re ->
    try
      if n = 0 then results,errors,n else (
      Printf.printf "%s %s\n%!" re.rxn_id re.reaction_smile;
      let time1 = Sys.time () in
      Types.time := time1;
      let r,messages = Smiles.prepare_record_for_matching re [] in
      let r = MatchingExec.translate_aromaticity r in
      let r = MatchingExec.set_unbreakable r in
      let messages,solutions,quality,no_broken_bonds = (*try*) AtomMapping.broken_bonds max_int 6 messages [r]
        (*with e -> (Printexc.to_string e) :: messages,[],0,0*) in
      let time2 = Sys.time () in
      (time2 -. time1,Xlist.size solutions,quality,no_broken_bonds,re) :: results, errors, n-1)
    with e -> results, (Printexc.to_string e,re) :: errors, n-1) in
  print_results4 res_filename results;
  print_errors err_filename errors;
  ()

(*let _ = test_broken_bonds 1000 "results/results4_200_OrganicSyntheses_Apr7_2015.txt" "results/errors_200_OrganicSyntheses_Apr7_2015.txt" (Import.load_organic_syntheses2 ())
let _ = test_broken_bonds 1000 "results/results4_przyklady_trudne.txt" "results/errors_przyklady_trudne.txt" (Import.load_przyklady_trudne ())
let _ = test_broken_bonds 1000 "results/results4_przyklady_trudne2.txt" "results/errors_przyklady_trudne2.txt" (Import.load_przyklady_trudne2 ())
let _ = test_broken_bonds 1000 "results/results4_Mapping_algorithm_train_set1.txt" "results/errors_Mapping_algorithm_train_set1.txt" (Import.load_paths ())
let _ = test_broken_bonds 1000 "results/results4_250_OrganicSyntheses.txt" "results/errors_250_OrganicSyntheses.txt" (Import.load_organic_syntheses ())
let _ = test_broken_bonds 1000 "results/results4_cukry_5_OHinSub.txt" "results/errors_cukry_5_OHinSub.txt" (Import.load_reactions Import.w_5_ohinsub_filename)
let _ = test_broken_bonds 1000 "results/results4_cukry_5_OHinProd.txt" "results/errors_cukry_5_OHinProd.txt" (Import.load_reactions Import.w_5_ohinprod_filename)
let _ = test_broken_bonds 1000 "results/results4_cukry_6_OHinProd.txt" "results/errors_cukry_6_OHinProd.txt" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = test_broken_bonds 1000 "results/results4_cukry_6_OHinSub.txt" "results/errors_cukry_6_OHinSub.txt" (Import.load_reactions Import.w_6_ohinsub_filename)*)
(* let _ = test_broken_bonds 1000 "results/results4_MetaCyc.txt" "results/errors_MetaCyc.txt" (Import.load_metacyc ()) *)
(*let _ = test_broken_bonds 10000 "results/results4_RT.txt" "results/errors_RT.txt" (Import.load_rt ())
let _ = test_broken_bonds 10000 "results/results4_rxn_db.txt" "results/errors_rxn_db.txt" (Import.load_rxn ())*)

(* let _ = test_broken_bonds 1000 "results/results4_pamieciozerne.txt" "results/errors_pamieciozerne.txt" (Xlist.map pamieciozerne (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)
(* let _ = MatchingExec.map_atoms2_distr 3 true "pamieciozerne_" "results/pamieciozerne_mapped.tex" (Xlist.map pamieciozerne (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let rec get_first n rev = function
    [] -> List.rev rev
  | x :: l -> if n = 0 then List.rev rev else get_first (n-1) (x :: rev) l

let rec drop_first n l =
  if n = 0 || l = [] then l else drop_first (n-1) (List.tl l)

(*let _ = MatchingExec.map_atoms2_distr2 4 true 1 5 "test_" "results/test_mapped_1-10.tex" (get_first 10 [] (Import.load_rt ()))
let _ = MatchingExec.map_atoms2_distr2 4 true 11 5 "test_" "results/test_mapped_11-20.tex" (get_first 10 [] (drop_first 10 (Import.load_rt ())))
let _ = MatchingExec.map_atoms2_distr2 4 true 21 5 "test_" "results/test_mapped_21-30.tex" (get_first 10 [] (drop_first 20 (Import.load_rt ())))*)

(* let rt = Import.load_rt ()
(*let _ = MatchingExec.map_atoms2_distr2 64 true 1 10000 "RT_" "results/RT_1_100000.tex" (get_first 100000 [] rt)*)
let _ = MatchingExec.map_atoms2_distr2 64 true 100001 10000 "RT_" "results/RT_100001_200000.tex" (get_first 100000 [] (drop_first 100000 rt))
let _ = MatchingExec.map_atoms2_distr2 64 true 200001 10000 "RT_" "results/RT_200001_300000.tex" (get_first 100000 [] (drop_first 200000 rt))
let _ = MatchingExec.map_atoms2_distr2 64 true 300001 10000 "RT_" "results/RT_300001_400000.tex" (get_first 100000 [] (drop_first 300000 rt))
let _ = MatchingExec.map_atoms2_distr2 64 true 400001 10000 "RT_" "results/RT_400001_500000.tex" (get_first 100000 [] (drop_first 400000 rt))
let _ = MatchingExec.map_atoms2_distr2 64 true 500001 10000 "RT_" "results/RT_500001_600000.tex" (get_first 100000 [] (drop_first 500000 rt))
let _ = MatchingExec.map_atoms2_distr2 64 true 600001 10000 "RT_" "results/RT_600001_700000.tex" (get_first 100000 [] (drop_first 600000 rt))
let rxn = Import.load_rxn ()
(* let _ = MatchingExec.map_atoms2_distr2 64 true 1 10000 "rxn_db_" "results/rxn_db_1_100000.tex" (get_first 100000 [] rxn) *)
let _ = MatchingExec.map_atoms2_distr2 64 true 100001 10000 "rxn_db_" "results/rxn_db_100001_200000.tex" (get_first 100000 [] (drop_first 100000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 200001 10000 "rxn_db_" "results/rxn_db_200001_300000.tex" (get_first 100000 [] (drop_first 200000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 300001 10000 "rxn_db_" "results/rxn_db_300001_400000.tex" (get_first 100000 [] (drop_first 300000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 400001 10000 "rxn_db_" "results/rxn_db_400001_500000.tex" (get_first 100000 [] (drop_first 400000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 500001 10000 "rxn_db_" "results/rxn_db_500001_600000.tex" (get_first 100000 [] (drop_first 500000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 600001 10000 "rxn_db_" "results/rxn_db_600001_700000.tex" (get_first 100000 [] (drop_first 600000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 700001 10000 "rxn_db_" "results/rxn_db_700001_800000.tex" (get_first 100000 [] (drop_first 700000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 800001 10000 "rxn_db_" "results/rxn_db_800001_900000.tex" (get_first 100000 [] (drop_first 800000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 900001 10000 "rxn_db_" "results/rxn_db_900001_1000000.tex" (get_first 100000 [] (drop_first 900000 rxn))
let _ = MatchingExec.map_atoms2_distr2 64 true 1000001 10000 "rxn_db_" "results/rxn_db_1000001_1100000.tex" (get_first 100000 [] (drop_first 1000000 rxn)) *)

(*
ssh wjaworski@students.mimuw.edu.pl
ssh wloczykij
cd Dropbox/Chemia/4efektywne_maczowanie/

http://bioputer.mimuw.edu.pl:8080/marvinjs/
*)

let pozostale_sierpien = [
  (* "36", "CS(=O)(Cl)=O.CCN(CC)CC.CC(C)([C@H](NC(C1=C(C=CC=C1)Br)=O)CO)C>>CC(C)([C@H]2COC(C3=C(C=CC=C3)Br)=N2)C.CS(=O)(O)=O.CC[NH+](CC)CC.[Cl-]";
  "121", "CC(OC(C)(C)C)=O.ClC1=C(C=CC(C#N)=C1)OC2=C(C=CC=C2)Br.O>>CC(C)(C)NC(C3=CC(Cl)=C(C=C3)OC4=C(C=CC=C4)Br)=O.CC(O)=O"; *)
  (* "126", "O=C(OC(C)=O)C.OC(/C=C\\C(NC1=CC=CC=C1)=O)=O>>O=C2N(C3=CC=CC=C3)C(C=C2)=O.O=C(O)C.O=C(O)C"; *)
  (* "179", "CCCCCCC(C=C)OC(CC(C)=O)=O>>CCCCCC/C=C/CC(C(C)=O)C(O)=O";*)
  "23", "C=CCOc1ccccc1>>Oc1ccccc1CC=C";
(*  "37", "CC#CC(O)=O.O=CC1=CC=CS1.CC(C)(C)[N+]#[C-].COc1ccc(CN)cc1>>COc1ccc(CN(C(C(=O)NC(C)(C)C)C2=CC=CS2)C(=O)C#CC)cc1.O";
  "38", "OC(C(F)(Sc1ccccc1)F)=O.O=Cc2ccccc2.[C-]#[N+]c3ccccc3.Nc4ccccc4>>FC(Sc5ccccc5)(C(N(c6ccccc6)C(c7ccccc7)C(Nc8ccccc8)=O)=O)F.O";
  "13", "COc1ccc(cc1)C#Cc1ccc(OC)cc1.C=C>>COc1ccc(cc1)C(=C)C(=C)c1ccc(OC)cc1";
  "40", "O=C1CCCCC1.CC(=O)\\C=C\\c1ccccc1>>O=C1CC(C2CCCCC2=C1)c1ccccc1"; *)
  ]

(* let _ = MatchingExec.print_mapped_reactions_rearr false "pozostale_sierpien_" "results/pozostale_sierpien.tex" (Xlist.map pozostale_sierpien (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let listopad = [
(* 200 Organic Syntheses *)
  "7", "CC(OCCCCCOC(C)=O)=O>>CC(O)=O.C=CCC=C.OC(C)=O";
  "94", "CC(OC(C)=O)=O.CC(OC(C)=O)=O.O=Cc1ccccc1.OC(CCC(c2ccccc2)=O)=O>>O=C3OC(c4ccccc4)=C/C3=C\\c5ccccc5.CC(O)=O.CC(O)=O.CC(O)=O.CC(O)=O";
(* Train set 1 *)
  "Path 3-5", "OC1=NC2=CC=CC=C2N=C1CC.CC(C)C=O.NCCOC.[C-]#[N+]C1CCCCC1>>CCC1=C(N(C(C(NC2CCCCC2)=O)C(C)C)CCOC)N=C3C(C=CC=C3)=N1.[H]O[H]";
  "Path 5-1", "CC(OC)(C)/C=C/C1=C[C@@H](O)[C@H](O2)[C@H]2C1=O.CC(OC)(C)/C=C/C1=C[C@@H](O)[C@H](O2)[C@H]2C1=O>>OC([C@H]1[C@@H]2O1)[C@@]3([H])[C@](C4=C[C@H](C(C)(OC)C)[C@@]3(/C=C/C(C)(OC)C)C2=O)([H])[C@@H](O)[C@H]5[C@H](O5)C4=O";
  "Path 8-13", "[H][C@]1([C@@H](CCCO)O)NC[C@H]2OC(C)(O[C@@H]12)C>>[H][C@@]12[C@@H]3OC(C)(O[C@@H]3CN1CCC[C@H]2O)C.O";
(* Przykłady trudne 1 *)
  "12", "CO[C@H]1C[C@H](S[Si](C)(C)C)[C@@H](CO1)OC(=O)N[Si](C)(C)C.CCc1cnc(O[Si](C)(C)C)nc1O[Si](C)(C)C>>CCC1=CN(C(CC2SC2CO)OC)C(=O)NC1=O.C[Si+](C)C.C[Si+](C)C.C[Si+](C)C.C[Si](C)(C)N.O=C=O";
  "17", "O=C1C(=O)c2cccc3cccc1c23.CC(=O)CC(=O)CC(C)=O.C1C2C=CC1C=C2>>CC(=O)c1ccc(C(C)=O)c2-c3cccc4cccc(-c12)c34.[C-]#[O+].C1C=CC=C1.O.O";
(* Przykłady trudne 2 *)
  "20", "CCOC(=O)CCCCC(=O)OCC>>O=C1CCCC1";
  "22", "Cc1ccc(cc1)S(=O)(=O)N(CC=C)C1CC(C=C1)N(CC=C)S(=O)(=O)c1ccc(C)cc1>>CC1=CC=C(C=C1)S(=O)(=O)N1CC=CC1C1C=CCN1S(=O)(=O)C1=CC=C(C)C=C1";
  ]

let ostatnia = [
  (* "0","[H][C@]12CC=C[C@@]1([H])CC(=C2)[C@@]1(CCC[C@H]1CC(OC)OC)O[Si](CC)(CC)CC>>[H][C@]12CC(OC)[C@]3([H])[C@@]4([H])CC=C[C@@]4([H])C[C@]13C(=O)CCC2.CC[Si](O)(CC)CC" *)
  "0","[H][C@]1([C@@H](CCCO)O)NC[C@H]2OC(C)(O[C@@H]12)C>>[H][C@@]12[C@@H]3OC(C)(O[C@@H]3CN1CCC[C@H]2O)C.O"
]

(* let _ = MatchingExec.print_mapped_reactions_rearr false "listopad_" "results/listopad.tex" (Xlist.map listopad (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)
(* let _ = MatchingExec.print_mapped_reactions_rearr false "ostatnia_" "results/ostatnia.tex" (Xlist.map ostatnia (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let macz_bands = [
  "1", "C(C1OC(OP([O-])([O-])=O)C(O)C1O)O.C12NC(=O)NC(=O)C=1NC(=O)N2>>C(C1C(O)C(O)C(N2C(=O)NC3C(NC(NC2=3)=O)=O)O1)O.[O-]P([O-])(O)=O";
  ]

(* let _ = MatchingExec.print_mapped_reactions_rearr false "macz_bands_" "results/macz_bands.tex" (Xlist.map macz_bands (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let estryfikacja = [
  "Path-2-1", "N[C@@H](CC1=CC=CC=C1)C(O)=O.CO>>N[C@@H](CC1=CC=CC=C1)C(OC)=O.[H]O[H]";
  "Path-2-1a", "N[C@@H](CC1=CC=CC=C1)C(O)=O>>N[C@@H](CC1=CC=CC=C1)C(OC)=O.[H]O[H]";
  "Path-2-1b", "N[C@@H](CC1=CC=CC=C1)C(O)=O.CO>>N[C@@H](CC1=CC=CC=C1)C(OC)=O";
  "Path-2-1c", "N[C@@H](CC1=CC=CC=C1)C(O)=O>>N[C@@H](CC1=CC=CC=C1)C(OC)=O";
  "20597-mod", "COC(=O)C=NNC1=Nc2ccc(Cl)cc2C(c2ccccc2)=NC1>>O=C1C=NN=C2N1-c1ccc(Cl)cc1C(c1ccccc1)=NC2";
  "20597-mod2", "COC(=O)C=NNC1=Nc2ccc(Cl)cc2C(c2ccccc2)=NC1>>O=C1C=NN=C2N1-c1ccc(Cl)cc1C(c1ccccc1)=NC2.CO";
  "20597-mod3", "COC(=O)C=NNC1=Nc2ccc(Cl)cc2C(c2ccccc2)=NC1>>CO";
  ] @ org_syn_estry

let patenty = [
(*  (* (Trans)Esterification *)
  "47494", "CCc1c(C=Cc2ccccc2)nc2sc(C(=O)O)cn2c1=O>>CCc1c(C=Cc2ccccc2)nc2sc(C(=O)OC)cn2c1=O";
  "85379", "CCOC(=O)C=Cc1ccc(NCCCCCCCCCCCCCC[Si](C)(C)C)cc1>>C[Si](C)(C)CCCCCCCCCCCCCCNc1ccc(C=CC(=O)O)cc1";
  (* (Trans)Esterification 2 *)
  "38448", "CCCCCCC(C)(C)c1cc(OCc2ccccc2)c2c(c1)OC(C)(C)CC2C(=O)O.O=C(Oc1ccc([N+](=O)[O-])cc1)C(F)(F)F>>CCCCCCC(C)(C)c1cc(OCc2ccccc2)c2c(c1)OC(C)(C)CC2C(=O)Oc1ccc([N+](=O)[O-])cc1";
  (* Aminocostam 2 *)
  "26353", "O=C(CNC(=O)OCc1ccccc1)ON1C(=O)CCC1=O.CCCCOC(OCCCC)C(N)CCCNC(=N)N>>CCCCOC(OCCCC)C(CCCNC(=N)N)NC(=O)CNC(=O)OCc1ccccc1";
  "58202", "CNCC(OC)OC.CC(C)c1cc(NC(=O)Oc2ccccc2)no1>>COC(CN(C)C(=O)Nc1cc(C(C)C)on1)OC";
  "69319", "CC1(C)OCC(CC=CCCCC(=O)O)C(c2ccccc2)O1>>CC1(C)OCC(CC=CCCCC(=O)NS(C)(=O)=O)C(c2ccccc2)O1";*)
  (* FIXME: find aux solvents *)
  (* FIXME: chlor *)
  (*"34665", "CCOC(=O)Cn1ncc(-c2ccccc2)c1-c1ccccc1.CN(C)CCCN>>Cl.CN(C)CCCNC(=O)Cn1ncc(-c2ccccc2)c1-c1ccccc1";
  "31131", "Cl.COC(=O)c1cc(C(O)CN2CCN(c3ccccc3OC)CC2)ccc1OC>>Cl.CNC(=O)c1cc(C(O)CN2CCN(c3ccccc3OC)CC2)ccc1OC";*)

  (* Ester reduction *)
  (*"9757", "CCOC(=O)c1ccc(COc2cccnc2N)cc1>>Nc1ncccc1OCc1ccc(CO)cc1";
  (*!*)"17206", "CCOC(=O)C(C)Oc1c([N+](=O)[O-])ccc(Oc2ccc([N+](=O)[O-])c(OC(C)C(=O)OCC)c2-c2ccc(C(F)(F)F)cc2Cl)c1-c1ccc(C(F)(F)F)cc1Cl>>CC(CO)Oc1c([N+](=O)[O-])ccc(Oc2ccc([N+](=O)[O-])c(OC(C)CO)c2-c2ccc(C(F)(F)F)cc2Cl)c1-c1ccc(C(F)(F)F)cc1Cl";*)
  (* Claisen Condensation *)
  (* "40584", "CC(=O)c1c(-c2ccccc2)cc2ccccn21.CCOC(=O)C(=O)OCC>>CCOC(=O)C(=O)CC(=O)c1c(-c2ccccc2)cc2ccccn21";
  "97584", "CCOC(=O)CC.COC(=O)c1ccc(OC)nn1>>CCOC(=O)C(C)C(=O)c1ccc(OC)nn1"; *)
  (* O-Alkilation of Enolates *)
  (* "44007", "CC(C)(C)C(=O)C(Oc1ccc(-c2ccc(Cl)cc2)cc1)n1ccnc1.CCOS(=O)(=O)OCC>>CCOC(=C(Oc1ccc(-c2ccc(Cl)cc2)cc1)n1ccnc1)C(C)(C)C"; *)
  (* Amidation 2 *)
  (* "50198", "CCOC(N)=O.Nc1ccc(Cl)cc1>>CCOC(=O)Nc1ccc(Cl)cc1"; *)
  (* Thioether Synthesis *)
  "56955", "COC1(NC(=O)CCCC(N)C(=O)O)C(=O)N2C(C(=O)O)=C(COC(C)=O)CSC21>>CCSCC1=C(C(=O)O)N2C(=O)C(NC(=O)CCCC(N)C(=O)O)(OC)C2SC1";
  (* Condensation with Phosphorus *)
  (* "88863", "COP(C)(=O)OC.CCCCC1(C(=O)OCC)CC1>>CCCCC1(C(=O)CP(=O)(OC)OC)CC1"; *)

  (*"20597", "COC(=O)C=NNC1=Nc2ccc(Cl)cc2C(c2ccccc2)=NC1>>O=c1cnnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2";
  "36827", "CC#N.c1cc2ccccc2s1.CCCCON=O>>N#CC(=NO)c1csc2ccccc12";
  "42059", "C=C(C)C(C(=O)OC)N1C(=O)C(N2C(=O)c3ccccc3C2=O)C1S(=O)Cl>>C=C1CS(=O)C2C(N3C(=O)c4ccccc4C3=O)C(=O)N2C1C(=O)OC";
  "71274", "CC(C#N)c1ccc(O)c(N)c1.O=C(O)C(O)c1ccccc1>>CC(C#N)c1ccc2oc(C(O)c3ccccc3)nc2c1";
  "81587", "CCOC(=O)C(=NOC)C(=O)CBr>>CCOC(=O)C(=NOC)c1coc(N)n1";
  "81931", "CN(N=O)C1=Nc2ccsc2C(c2ccccc2)=NC1>>O=[N+]([O-])C=C1CN=C(c2ccccc2)c2sccc2N1";
  (*błąd w danych*)"90520", "CS(=O)(=O)C(O)C1COC(c2ccccc2)=N1.CCCCCCCCCCCCCCCCCCS>>CCCCCCCCCCCCCCCCCCSCC1COC(c2ccccc2)=N1";
  "97699", "CSC(=C[N+](=O)[O-])S(C)=O.Cc1ncnc1CSCCN>>Cc1ncnc1CSCCNC(=C[N+](=O)[O-])NCC(F)(F)F";*)
  ]

let bledy = [
  "38448", "CCCCCCC(C)(C)c1cc(OCc2ccccc2)c2c(c1)OC(C)(C)CC2C(=O)O.O=C(Oc1ccc([N+](=O)[O-])cc1)C(F)(F)F>>CCCCCCC(C)(C)c1cc(OCc2ccccc2)c2c(c1)OC(C)(C)CC2C(=O)Oc1ccc([N+](=O)[O-])cc1";
  "Path-2-3", "O=C(OC)[C@@H](NCC=C)CC1=CC=CC=C1.C#CCBr>>O=C(OC)[C@@H](N(CC#C)CC=C)CC1=CC=CC=C1.[H]Br";
  (* "Path-2-4", "O=C(OC)[C@@H](N(CC#C)CC=C)CC1=CC=CC=C1>>O=C(OC)[C@@H](N1CC(C=C)=CC1)CC2=CC=CC=C2"; *)
(*  "Path-6-2", "O[C@@H](CC=C)C1=CC=CC=C1.ClC(N(C(C)C)C(C)C)=O>>C=CC[C@H](OC(N(C(C)C)C(C)C)=O)C1=CC=CC=C1.[H]Cl";
  "Path-7-1", "CC1=C2[C@@](CCC2)([H])[C@]3([H])[C@@](/C(C(O3)=O)=C/C)([H])[C@@H](O)C1>>O[C@H]1CC(C)=C2[C@@]([C@]1([H])/C3=C/C)([H])[C@@](CCC2)([H])OC3=O";
  "Path-8-3", "CC1(O[C@@H]2COC([C@@H]2O1)=O)C>>CC1(O[C@@H]2COC([C@@H]2O1)O)C";
  "Path-8-8", "COC(CC/C=C/[C@@H]1OC(C)(O[C@@H]1CO)C)=O>>COC(CC/C=C/[C@@H]1OC(C)(O[C@@H]1C=O)C)=O";
  "Path-8-9", "COC(CC/C=C/[C@@H]1OC(C)(O[C@@H]1C=O)C)=O.O>>COC(CC/C=C/[C@@H]1OC(C)(O[C@@H]1C(O)=O)C)=O";*)
  ]

let ritter = [
  "31", "C=CC#N.OCC1=CC=CC=C1>>C=CC(=O)NCC1=CC=CC=C1";
  "109", "CC(OC(C)(C)C)=O.COC(C1=CC=C(C#N)C=C1)=O.O>>COC(C2=CC=C(C(NC(C)(C)C)=O)C=C2)=O.CC(O)=O";
  "121", "CC(OC(C)(C)C)=O.ClC1=C(C=CC(C#N)=C1)OC2=C(C=CC=C2)Br.O>>CC(C)(C)NC(C3=CC(Cl)=C(C=C3)OC4=C(C=CC=C4)Br)=O.CC(O)=O";
  "122", "CC(OC(C)(C)C)=O.N#CCCC1=CC=CC=C1.O>>CC(O)=O.CC(C)(NC(CCC2=CC=CC=C2)=O)C";
  "123", "CC(OC(C)(C)C)=O.CC(NC1=CC=C(C#N)C=C1)=O.O>>CC(NC2=CC=C(C(NC(C)(C)C)=O)C=C2)=O.CC(O)=O";
  ]

(* let _ = MatchingExec.print_mapped_reactions_rearr2 "estyfikacja_" "results/estryfikacja.tex" (Xlist.map estryfikacja (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)
(* let _ = MatchingExec.print_mapped_reactions_rearr2 "perycykliczne_" "results/perycykliczne.tex"
   [{empty_record with reaction_smile = "C\\C=C\\[C@H](O)C\\C=C\\c1ccccc1>>C[C@H](CC=O)[C@@H](C=C)c1ccccc1"}] *)
(* let _ = MatchingExec.print_mapped_reactions_rearr2 "patenty_" "results/patenty.tex" (Xlist.map patenty (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)
(* let _ = MatchingExec.print_mapped_reactions_rearr2 "ritter_" "results/ritter.tex" (Xlist.map ritter (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)
(* let _ = MatchingExec.print_mapped_reactions_rearr2 "bledy_" "results/bledy.tex" (Xlist.map bledy (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let ostatnie = [
  "1", "CC(OC)(C)/C=C/C1=C[C@@H](O)[C@H](O2)[C@H]2C1=O.CC(OC)(C)/C=C/C1=C[C@@H](O)[C@H](O2)[C@H]2C1=O>>OC([C@H]1[C@@H]2O1)[C@@]3([H])[C@](C4=C[C@H](C(C)(OC)C)[C@@]3(/C=C/C(C)(OC)C)C2=O)([H])[C@@H](O)[C@H]5[C@H](O5)C4=O";
  "2", "C=C/C=C/OC(C)=O.C1=CC=CCC=C1>>[H][C@@]23C=CC=C[C@@](CC=C[C@@H]3OC(C)=O)([H])C2";
  "3", "CC(OC(C)=O)=O.CC(OC(C)=O)=O.O=Cc1ccccc1.OC(CCC(c2ccccc2)=O)=O>>O=C3OC(c4ccccc4)=C/C3=C\\c5ccccc5.CC(O)=O.CC(O)=O.CC(O)=O.CC(O)=O";
  "4", "C1N=C1c1ccccc1.O=C1C(=C(C(=C1c1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1>>c1ccc(cc1)C1C(=NC(c2ccccc2)=C(c2ccccc2)C(c2ccccc2)=C1c1ccccc1)c1ccccc1";
  "5", "[H][C@]12CC=C[C@@]1([H])CC(=C2)[C@@]1(CCC[C@H]1CC(OC)OC)O[Si](CC)(CC)CC>>[H][C@]12CC(OC)[C@]3([H])[C@@]4([H])CC=C[C@@]4([H])C[C@]13C(=O)CCC2.CC[Si](O)(CC)CC";
  "6", "Cc1ccc(cc1)S(=O)(=O)N(CC=C)C1CC(C=C1)N(CC=C)S(=O)(=O)c1ccc(C)cc1>>CC1=CC=C(C=C1)S(=O)(=O)N1CC=CC1C1C=CCN1S(=O)(=O)C1=CC=C(C)C=C1";
  "7", "O=C1C(=O)c2cccc3cccc1c23.CC(=O)CC(=O)CC(C)=O.C1C2C=CC1C=C2>>CC(=O)c1ccc(C(C)=O)c2-c3cccc4cccc(-c12)c34.[C-]#[O+].C1C=CC=C1.O.O";
  "8", "CO[C@H]1C[C@H](S[Si](C)(C)C)[C@@H](CO1)OC(=O)N[Si](C)(C)C.CCc1cnc(O[Si](C)(C)C)nc1O[Si](C)(C)C>>CCC1=CN(C(CC2SC2CO)OC)C(=O)NC1=O.C[Si+](C)C.C[Si+](C)C.C[Si+](C)C.C[Si](C)(C)N.O=C=O";
  "9", "OC1=NC2=CC=CC=C2N=C1CC.CC(C)C=O.NCCOC.[C-]#[N+]C1CCCCC1>>CCC1=C(N(C(C(NC2CCCCC2)=O)C(C)C)CCOC)N=C3C(C=CC=C3)=N1.[H]O[H]";
  "10", "[H][C@]1([C@@H](CCCO)O)NC[C@H]2OC(C)(O[C@@H]12)C>>[H][C@@]12[C@@H]3OC(C)(O[C@@H]3CN1CCC[C@H]2O)C.O";
  "11", "NC(CO)C1=CC=CC=C1.ClS(Cl)=O>>[NH3+]C(CCl)C2=CC=CC=C2.O=S=O.[Cl-]";
  "12", "ClS(Cl)=O.NCCC1=C(C=C(C=C1)Cl)CO.[OH-].[OH-]>>ClC2=CC=C3CCNCC3=C2.O=S=O.O.O.[Cl-].[Cl-]";
  "13", "CC(OCCCCCOC(C)=O)=O>>CC(O)=O.C=CCC=C.OC(C)=O";
]

(* let _ = MatchingExec.print_mapped_reactions_rearr2 "losowe_patenty_" "results/losowe_patenty.tex" (Import.load_losowe_patenty ()) *)
(* let _ = MatchingExec.print_mapped_reactions_rearr2 "train_set1_" "results/Mapping_algorithm_train_set1.tex" (Import.load_paths ()) *)
(*let _ = MatchingExec.print_mapped_reactions_rearr2 "200_OrganicSyntheses_Apr7_2015_" "results/200_OrganicSyntheses_Apr7_2015.tex" (Import.load_organic_syntheses2 ()) *)
(*let _ = MatchingExec.print_mapped_reactions_rearr2 "przyklady_trudne_" "results/przyklady_trudne.tex" (Import.load_przyklady_trudne ())
let _ = MatchingExec.print_mapped_reactions_rearr2 "przyklady_trudne2_" "results/przyklady_trudne2.tex" (Import.load_przyklady_trudne2 ())
let _ = MatchingExec.print_mapped_reactions_rearr2 "cukry_5_OHinProd_" "results/cukry_5_OHinProd.tex" (Import.load_reactions Import.w_5_ohinprod_filename)
let _ = MatchingExec.print_mapped_reactions_rearr2 "cukry_5_OHinSub_" "results/cukry_5_OHinSub.tex" (Import.load_reactions Import.w_5_ohinsub_filename)
let _ = MatchingExec.print_mapped_reactions_rearr2 "cukry_6_OHinProd_" "results/cukry_6_OHinProd.tex" (Import.load_reactions Import.w_6_ohinprod_filename)
let _ = MatchingExec.print_mapped_reactions_rearr2 "cukry_6_OHinSub_" "results/cukry_6_OHinSub.tex" (Import.load_reactions Import.w_6_ohinsub_filename)*)
(* let _ = MatchingExec.print_mapped_reactions_rearr false "ostatnie_" "results/ostatnie.tex" (Xlist.map ostatnie (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let sierpien2017 = [
  "1", "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\COP([O-])(=O)OP([O-])([O-])=O.CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\COP([O-])(=O)OP([O-])([O-])=O>>CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C.OP([O-])(=O)OP([O-])([O-])=O.OP([O-])(=O)OP([O-])([O-])=O";
  "2", "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C=C(/C)CC\\C=C(/C)CCC1OC1(C)C>>C[C@H](CCC=C(C)C)[C@H]1CC[C@]2(C)C1CCC1=C2CC[C@H]2C(C)(C)[C@@H](O)CC[C@]12C";
]

(* let _ = MatchingExec.print_mapped_reactions_rearr false "sierpien2017_" "results/sierpien2017.tex" (Xlist.map sierpien2017 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let listopad2017 = [
  "1", "CO[C@H]1C[C@H](S[Si](C)(C)C)[C@@H](CO1)OC(=O)N[Si](C)(C)C.CCc1cnc(O[Si](C)(C)C)nc1O[Si](C)(C)C>>CCC1=CN(C(CC2SC2CO)OC)C(=O)NC1=O.C[Si+](C)C.C[Si+](C)C.C[Si+](C)C.C[Si](C)(C)N.O=C=O";
  ]

(* let _ = MatchingExec.print_mapped_reactions_rearr false "listopad2017_" "results/listopad2017.tex" (Xlist.map listopad2017 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)

let kwiecien2018 = [
  "A14", "C1=CC=C(C=C1)P(C1=CC=CC=C1)C1=CC=CC=C1.COC1=CC=C(COC\\C=C(/C)[C@@H]2OC3(CCCC3)O[C@H]2CO)C=C1.II>>COC1=CC=C(COC\\C=C(/C)[C@@H]2OC3(CCCC3)O[C@H]2CI)C=C1.O=P(C1=CC=CC=C1)(C1=CC=CC=C1)C1=CC=CC=C1";
  "A148", "OCCc1ccc2OCCc2c1>>BrCCc1ccc2OCCc2c1";
  "A176", "COc1ccc(N)c(N)c1.OC(=O)C(O)=O>>COc1ccc2nc(O)c(O)nc2c1";
  "112", "CC1=CC=C(C=C1)C(=O)C(O)O.CC1=CC=C(C=C1)C(=O)C(O)O.NCCC(O)=O.CC1(C)CC(=O)C=C(C1)NC1=CC=C(F)C=C1>>[H][C@]1(OC2(N(C3=C(C(=O)CC(C)(C)C3)[C@]2([H])N1CCC(O)=O)C1=CC=C(F)C=C1)C1=CC=C(C)C=C1)C(=O)C1=CC=C(C)C=C1";
  "117", "COC(=O)C(=C(/[S-])C1=CC=C(C=C1)N(C)C)\\C1=CC=C(OC)C=C1.CCOC(=O)C1=C(CBr)C=CC=C1>>CCOC(=O)C1=CC=CC=C1C1=C(O)C(=C(S1)C1=CC=C(C=C1)N(C)C)C1=CC=C(OC)C=C1";
  "120", "CC(C)(C)[Si](C)(C)O[C@@H](CC[C@@H](CC=C)OS(C)(=O)=O)COCC1=CC=CC=C1>>C=CC[C@H]1CC[C@@H](COCC2=CC=CC=C2)O1";
  "136", "BrC1=CC=C(\\C=C\\C(=O)N2C=CC=N2)C=C1.CCOC(=O)C\\N=C\\C1=CC=CC=C1.O=C1C=CC(=O)N1C1=CC=CC=C1.O=C\\C=C\\C1=CC=CC=C1>>CCOC(=O)[C@]12[C@H]3[C@@H]([C@@H](\\C=C\\C4=CC=CC=C4)N1[C@@H]([C@@H]([C@H]2C1=CC=C(Br)C=C1)C(=O)N1C=CC=N1)C1=CC=CC=C1)C(=O)N(C3=O)C1=CC=CC=C1";
  "137", "COC1=C(Br)C(C(OC(=O)C2=CCC2)C=C)=C(OC)C=C1>>COC1=C(Br)C(C2OC(=O)C(CCC=C)=C2)=C(OC)C=C1";
  "143", "CC(C)CC(=O)C1=C(O)C=C(O)C=C1O.COC(CC(C)C)OC.CC1(C)C(O)=CC(=O)C(C)(C)C1=O>>CC(C)CC1C2=C(OC3=C1C(=O)C(C)(C)C(=O)C3(C)C)C(C(=O)CC(C)C)=C(O)C=C2O";
  "147", "BrC1=CC=CC2=CC=CC=C12.C1=CC=C(C=C1)C#CC1=CC=CC=C1>>C1C2=CC=CC=C2C2=C1C1=CC=CC=C1C=C2";
  "161", "CC(C)(C)OC(=O)N(CC=C)[C@@H]1CC=C[C@H]1CN(CC=C)S(=O)(=O)C1=CC=C(C=C1)[N+]([O-])=O>>[H][C@]1(CN(CC=C1)S(=O)(=O)C1=CC=C(C=C1)[N+]([O-])=O)[C@@]1([H])CC=CCN1C(=O)OC(C)(C)C";
  "167", "CC[C@@H](O[Si](C)(C)C(C)(C)C)[C@@H](O[C@@H](CCOCOCCOC)C(=O)CS(=O)(=O)C1=CC=CC=C1)\\C=C\\COC(C)=O>>CC[C@@H](O[Si](C)(C)C(C)(C)C)[C@H]1O[C@@H](CCOCOCCOC)C(=O)CC\\C=C\\1";
  "168", "COCOC[C@@]1(C)[C@@H]2C[C@@H](C=C2)[C@H]1C(=O)C(C)=C.C=C>>[H][C@@]12C[C@H](C=C)[C@H](COCOC)[C@]1([H])C(=O)C(C)=C2";
  "172", "C\\C=C1\\CN2CCC3=C(NC4=C3C=CC=C4)[C@]2(CO)[C@H]1CCO>>[H][C@]12CCN(C\\C1=C\\C)CCC1=C(NC3=C1C=CC=C3)C2=C"
  ]

(*let _ = MatchingExec.print_mapped_reactions_rearr false "kwiecien2018_" "results/kwiecien2018.tex" (Xlist.map kwiecien2018 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s}))*)

(*let _ =
  let l = Import.load_wszystkie () in
  Xlist.iter l (fun r -> Printf.printf "%d %s\n" (Smiles.count_reactants_atoms r) r.reaction_smile)*)

let listopad2018 = [
  "1", "OCCCC1=CC=CC=C1.BrP(Br)Br>>BrCCCC1=CC=CC=C1.OP(Br)Br";
(*  "1a", "OCCCC.BrP(Br)Br>>BrCCCC.OP(Br)Br";*)
(*   "1b", "OCCCC.BrP>>BrCCCC.OP"; *)
(*    "2", "CC1(C)C2CCC1(C)C(O)C2>>CC1(C)C2CCC(C2)C1=C.O"; *)
(*  "2", "O=C1OC(=O)C2=CC3=C(C=C12)C(=O)OC3=O.C1=CC=CC=C1.C1=CC=CC=C1>>OC(=O)C1=CC(C(=O)C2=CC=CC=C2)=C(C=C1C(=O)C1=CC=CC=C1)C(O)=O";*)
  "3", "OCC1=CC2=CC=CC=C2C=C1.BrS(Br)=O>>OS(Br)=O.BrCC1=CC2=CC=CC=C2C=C1";
  "4", "O=[Se]=O>>CN(C)S(=O)(=O)C1=C2C=CC(C)=NC2=C(O)C=C1.CN(C)S(=O)(=O)C1=C2C=CC(C=O)=NC2=C(O)C=C1.O=[Se]";
  "5", "CCCCCCO.[Br-]>>CCCCCCBr.[OH-]";
  "6", "CCC(=O)OC.I>>CCC(O)=O.CI";
  "7", "CC(=O)OCCS>>CC(O)=O.C1CS1";
  "8", "CC(C)[C@H](N)CO>>CC(C)[C@H]1CN1"
  ]

(* let _ = MatchingExec.print_labeled_reactions "listopad2018" "results/listopad2018.tex" (Xlist.map listopad2018 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)
(* let _ = MatchingExec.print_mapped_reactions "listopad2018" "results/listopad2018.tex" (Xlist.map listopad2018 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s})) *)
(*let _ = MatchingExec.print_mapped_reactions_rearr false "listopad2018" "results/listopad2018.tex" (Xlist.map listopad2018 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s}))*)
let _ = MatchingExec.print_mapped_reactions_rearr true "listopad2018" "results/listopad2018.tex" (Xlist.map listopad2018 (fun (id,s) -> {empty_record with rxn_id=id; reaction_smile=s}))
