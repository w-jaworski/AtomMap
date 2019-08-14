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
open Xstd

(* reprezentacje cząsteczki:
- smiles jako napis
- smiles jako drzewo (dodane etykiety wierzchołków) -> konwersja na napis
- drzewo z kątami i dodanymi wodorami -> wizualizacja
- graf
*)

type symbol =
    A of smile_tree
  | LP
  | RP
  | LB
  | RB
  | Number of string
  | Symb of string

(***********************************************************************************************)
(* diagnostyczna konwersja na string *)

let diagnostic_string_of_valence = function
    Alifatic -> "Alifatic"
  | Aromatic -> "Aromatic"

let diagnostic_string_of_chirality = function
    AntiClockwise -> "AntiClockwise"
  | Clockwise -> "Clockwise"
  | Unspecified -> "Unspecified"

let diagnostic_string_of_bond = function
    Default -> "Default"
  | Single -> "Single"
  | Slash -> "Slash"
  | Backslash -> "Backslash"
  | Double -> "Double"
  | Triple -> "Triple"
  | Aromatic -> "Aromatic"

let rec diagnostic_string_of_smile_tree = function
    SAtom(p,q) ->
      "Atom(" ^ p.sname ^ "," ^ diagnostic_string_of_valence p.svalence ^ "," ^ string_of_int p.scharge ^ "," ^
        diagnostic_string_of_chirality p.schirality ^ "," ^ string_of_int p.shydrogens ^ ",[" ^ String.concat "; " (Xlist.map q (fun (b,a) ->
        diagnostic_string_of_bond b ^ "," ^ diagnostic_string_of_smile_tree a)) ^ "]," ^ string_of_int p.sid ^ ")"
  | SLink s -> "Link(" ^ s ^ ")"

let diagnostic_string_of_smile_tree_simple = function
    SAtom(p,q) ->
      "Atom(" ^ p.sname ^ ")"
  | SLink s -> "Link(" ^ s ^ ")"

let rec diagnostic_string_of_symbol = function
    A a -> diagnostic_string_of_smile_tree a
  | LP -> "LP"
  | RP -> "RP"
  | LB -> "LB"
  | RB -> "RB"
  | Number s -> "Number " ^ s
  | Symb s -> "Symb " ^ s

let rec diagnostic_string_of_symbol_simple = function
    A a -> diagnostic_string_of_smile_tree_simple a
  | LP -> "LP"
  | RP -> "RP"
  | LB -> "LB"
  | RB -> "RB"
  | Number s -> "Number " ^ s
  | Symb s -> "Symb " ^ s

(***********************************************************************************************)
(* translacja ze smilesu jako napisu do smilesu jako drzewa *)

let make_atom name valence normality =
  SAtom({sname=name; svalence=valence; scharge=0; schirality=Unspecified; shydrogens=normality; sid=0},[])

let rec merge_atom_names id s = function
    "C" :: "l" :: l -> (A (make_atom "Cl" Alifatic (-1))) :: (merge_atom_names id s l)
  | "B" :: "r" :: l -> (A (make_atom "Br" Alifatic (-1))) :: (merge_atom_names id s l)
  | "[" :: "H" :: "e" :: l -> LB :: (A (make_atom "He" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "L" :: "i" :: l -> LB :: (A (make_atom "Li" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "B" :: "e" :: l -> LB :: (A (make_atom "Be" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: "e" :: l -> LB :: (A (make_atom "Ne" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: "a" :: l -> LB :: (A (make_atom "Na" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "M" :: "g" :: l -> LB :: (A (make_atom "Mg" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "l" :: l -> LB :: (A (make_atom "Al" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: "i" :: l -> LB :: (A (make_atom "Si" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "l" :: l -> LB :: (A (make_atom "Cl" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "r" :: l -> LB :: (A (make_atom "Ar" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "a" :: l -> LB :: (A (make_atom "Ca" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: "c" :: l -> LB :: (A (make_atom "Sc" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "i" :: l -> LB :: (A (make_atom "Ti" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "r" :: l -> LB :: (A (make_atom "Cr" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "M" :: "n" :: l -> LB :: (A (make_atom "Mn" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "F" :: "e" :: l -> LB :: (A (make_atom "Fe" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "o" :: l -> LB :: (A (make_atom "Co" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: "i" :: l -> LB :: (A (make_atom "Ni" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "u" :: l -> LB :: (A (make_atom "Cu" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "Z" :: "n" :: l -> LB :: (A (make_atom "Zn" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "G" :: "a" :: l -> LB :: (A (make_atom "Ga" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "G" :: "e" :: l -> LB :: (A (make_atom "Ge" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "s" :: l -> LB :: (A (make_atom "As" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: "e" :: l -> LB :: (A (make_atom "Se" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "B" :: "r" :: l -> LB :: (A (make_atom "Br" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "K" :: "r" :: l -> LB :: (A (make_atom "Kr" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "b" :: l -> LB :: (A (make_atom "Rb" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: "r" :: l -> LB :: (A (make_atom "Sr" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "Z" :: "r" :: l -> LB :: (A (make_atom "Zr" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: "b" :: l -> LB :: (A (make_atom "Nb" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "M" :: "o" :: l -> LB :: (A (make_atom "Mo" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "c" :: l -> LB :: (A (make_atom "Tc" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "u" :: l -> LB :: (A (make_atom "Ru" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "h" :: l -> LB :: (A (make_atom "Rh" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "d" :: l -> LB :: (A (make_atom "Pd" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "g" :: l -> LB :: (A (make_atom "Ag" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "d" :: l -> LB :: (A (make_atom "Cd" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "I" :: "n" :: l -> LB :: (A (make_atom "In" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: "n" :: l -> LB :: (A (make_atom "Sn" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: "b" :: l -> LB :: (A (make_atom "Sb" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "e" :: l -> LB :: (A (make_atom "Te" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "X" :: "e" :: l -> LB :: (A (make_atom "Xe" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "s" :: l -> LB :: (A (make_atom "Cs" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "B" :: "a" :: l -> LB :: (A (make_atom "Ba" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "H" :: "f" :: l -> LB :: (A (make_atom "Hf" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "a" :: l -> LB :: (A (make_atom "Ta" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "e" :: l -> LB :: (A (make_atom "Re" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "O" :: "s" :: l -> LB :: (A (make_atom "Os" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "I" :: "r" :: l -> LB :: (A (make_atom "Ir" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "t" :: l -> LB :: (A (make_atom "Pt" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "u" :: l -> LB :: (A (make_atom "Au" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "H" :: "g" :: l -> LB :: (A (make_atom "Hg" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "l" :: l -> LB :: (A (make_atom "Tl" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "b" :: l -> LB :: (A (make_atom "Pb" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "B" :: "i" :: l -> LB :: (A (make_atom "Bi" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "o" :: l -> LB :: (A (make_atom "Po" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "t" :: l -> LB :: (A (make_atom "At" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "n" :: l -> LB :: (A (make_atom "Rn" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "F" :: "r" :: l -> LB :: (A (make_atom "Fr" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "a" :: l -> LB :: (A (make_atom "Ra" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "L" :: "a" :: l -> LB :: (A (make_atom "La" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "e" :: l -> LB :: (A (make_atom "Ce" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "r" :: l -> LB :: (A (make_atom "Pr" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: "d" :: l -> LB :: (A (make_atom "Nd" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "m" :: l -> LB :: (A (make_atom "Pm" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: "m" :: l -> LB :: (A (make_atom "Sm" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "E" :: "u" :: l -> LB :: (A (make_atom "Eu" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "G" :: "d" :: l -> LB :: (A (make_atom "Gd" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "b" :: l -> LB :: (A (make_atom "Tb" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "D" :: "y" :: l -> LB :: (A (make_atom "Dy" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "H" :: "o" :: l -> LB :: (A (make_atom "Ho" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "E" :: "r" :: l -> LB :: (A (make_atom "Er" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "m" :: l -> LB :: (A (make_atom "Tm" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "Y" :: "b" :: l -> LB :: (A (make_atom "Yb" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "L" :: "u" :: l -> LB :: (A (make_atom "Lu" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "c" :: l -> LB :: (A (make_atom "Ac" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "T" :: "h" :: l -> LB :: (A (make_atom "Th" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "a" :: l -> LB :: (A (make_atom "Pa" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: "p" :: l -> LB :: (A (make_atom "Np" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: "u" :: l -> LB :: (A (make_atom "Pu" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "A" :: "m" :: l -> LB :: (A (make_atom "Am" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "m" :: l -> LB :: (A (make_atom "Cm" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "B" :: "k" :: l -> LB :: (A (make_atom "Bk" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: "f" :: l -> LB :: (A (make_atom "Cf" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "E" :: "s" :: l -> LB :: (A (make_atom "Es" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "F" :: "m" :: l -> LB :: (A (make_atom "Fm" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "M" :: "d" :: l -> LB :: (A (make_atom "Md" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: "o" :: l -> LB :: (A (make_atom "No" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "L" :: "r" :: l -> LB :: (A (make_atom "Lr" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "s" :: "e" :: l -> LB :: (A (make_atom "Se" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "t" :: "e" :: l -> LB :: (A (make_atom "Te" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "s" :: "i" :: l -> LB :: (A (make_atom "Si" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "H" :: l -> LB :: (A (make_atom "H" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "B" :: l -> LB :: (A (make_atom "B" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "C" :: l -> LB :: (A (make_atom "C" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "N" :: l -> LB :: (A (make_atom "N" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "O" :: l -> LB :: (A (make_atom "O" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "F" :: l -> LB :: (A (make_atom "F" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "P" :: l -> LB :: (A (make_atom "P" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "S" :: l -> LB :: (A (make_atom "S" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "K" :: l -> LB :: (A (make_atom "K" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "V" :: l -> LB :: (A (make_atom "V" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "Y" :: l -> LB :: (A (make_atom "Y" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "I" :: l -> LB :: (A (make_atom "I" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "W" :: l -> LB :: (A (make_atom "W" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "U" :: l -> LB :: (A (make_atom "U" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "b" :: l -> LB :: (A (make_atom "B" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "c" :: l -> LB :: (A (make_atom "C" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "n" :: l -> LB :: (A (make_atom "N" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "o" :: l -> LB :: (A (make_atom "O" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "s" :: l -> LB :: (A (make_atom "S" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "p" :: l -> LB :: (A (make_atom "P" Aromatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "1" :: l -> LB :: (A (make_atom "R1" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "2" :: l -> LB :: (A (make_atom "R2" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "3" :: l -> LB :: (A (make_atom "R3" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "4" :: l -> LB :: (A (make_atom "R4" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "5" :: l -> LB :: (A (make_atom "R5" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "6" :: l -> LB :: (A (make_atom "R6" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "7" :: l -> LB :: (A (make_atom "R7" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "8" :: l -> LB :: (A (make_atom "R8" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: "9" :: l -> LB :: (A (make_atom "R9" Alifatic 0)) :: (merge_atom_names id s l)
  | "[" :: "R" :: l -> LB :: (A (make_atom "R" Alifatic 0)) :: (merge_atom_names id s l)
  | "@" :: "@" :: "@" :: l -> failwith "merge_atom_names: ni"
  | "@" :: "@" :: l -> Symb "@@" :: (merge_atom_names id s l)
  | "B" :: l -> (A (make_atom "B" Alifatic (-1))) :: (merge_atom_names id s l)
  | "C" :: l -> (A (make_atom "C" Alifatic (-1))) :: (merge_atom_names id s l)
  | "N" :: l -> (A (make_atom "N" Alifatic (-1))) :: (merge_atom_names id s l)
  | "O" :: l -> (A (make_atom "O" Alifatic (-1))) :: (merge_atom_names id s l)
  | "P" :: l -> (A (make_atom "P" Alifatic (-1))) :: (merge_atom_names id s l)
  | "S" :: l -> (A (make_atom "S" Alifatic (-1))) :: (merge_atom_names id s l)
  | "H" :: l -> (A (make_atom "H" Alifatic (-1))) :: (merge_atom_names id s l)
  | "F" :: l -> (A (make_atom "F" Alifatic (-1))) :: (merge_atom_names id s l)
  | "I" :: l -> (A (make_atom "I" Alifatic (-1))) :: (merge_atom_names id s l)
  | "0" :: l -> Number "0" :: (merge_atom_names id s l)
  | "1" :: l -> Number "1" :: (merge_atom_names id s l)
  | "2" :: l -> Number "2" :: (merge_atom_names id s l)
  | "3" :: l -> Number "3" :: (merge_atom_names id s l)
  | "4" :: l -> Number "4" :: (merge_atom_names id s l)
  | "5" :: l -> Number "5" :: (merge_atom_names id s l)
  | "6" :: l -> Number "6" :: (merge_atom_names id s l)
  | "7" :: l -> Number "7" :: (merge_atom_names id s l)
  | "8" :: l -> Number "8" :: (merge_atom_names id s l)
  | "9" :: l -> Number "9" :: (merge_atom_names id s l)
  | "-" :: l -> Symb "-" :: (merge_atom_names id s l)
  | "+" :: l -> Symb "+" :: (merge_atom_names id s l)
  | "=" :: l -> Symb "=" :: (merge_atom_names id s l)
  | ":" :: l -> Symb ":" :: (merge_atom_names id s l)
  | "%" :: l -> Symb "%" :: (merge_atom_names id s l)
(*   | "[" :: l -> LB :: (merge_atom_names id s l) *)
  | "]" :: l -> RB :: (merge_atom_names id s l)
  | "(" :: l -> LP :: (merge_atom_names id s l)
  | ")" :: l -> RP :: (merge_atom_names id s l)
  | "b" :: l -> (A (make_atom "B" Aromatic (-1))) :: (merge_atom_names id s l)
  | "c" :: l -> (A (make_atom "C" Aromatic (-1))) :: (merge_atom_names id s l)
  | "n" :: l -> (A (make_atom "N" Aromatic (-1))) :: (merge_atom_names id s l)
  | "o" :: l -> (A (make_atom "O" Aromatic (-1))) :: (merge_atom_names id s l)
  | "s" :: l -> (A (make_atom "S" Aromatic (-1))) :: (merge_atom_names id s l)
  | "p" :: l -> (A (make_atom "P" Aromatic (-1))) :: (merge_atom_names id s l)
  | "#" :: l -> Symb "#" :: (merge_atom_names id s l)
  | "@" :: l -> Symb "@" :: (merge_atom_names id s l)
  | "/" :: l -> Symb "/" :: (merge_atom_names id s l)
  | "\\" :: l -> Symb "\\" :: (merge_atom_names id s l)
  | [] -> []
  | l -> failwith ("merge_atom_names " ^ id ^ " " ^ s ^ ": " ^ (String.concat "; " l))

let rec merge_numbers id s = function
    [] -> []
  | Symb "%" :: Number a :: Number b :: Number c :: Number d :: _ -> failwith "merge_numbers: ni"
  | Symb "%" :: Number a :: Number b :: Number c :: l -> Number(a ^ b ^ c) :: (merge_numbers id s l)
  | Symb "%" :: Number a :: Number b :: l -> Number(a ^ b) :: (merge_numbers id s l)
  | Symb "%" :: _ -> failwith "merge_numbers 1"
  | Symb ":" :: Number a :: Number b :: Number c :: Number d :: _ -> failwith "merge_numbers: ni"
  | Symb ":" :: Number a :: Number b :: Number c :: l -> Symb ":" :: Number(a ^ b ^ c) :: (merge_numbers id s l)
  | Symb ":" :: Number a :: Number b :: l -> Symb ":" :: Number(a ^ b) :: (merge_numbers id s l)
  | x :: l -> x :: (merge_numbers id s l)


let rec merge_atoms_inside id s = function
    SAtom(p,q), RB :: l -> SAtom(p,q),l
  | SAtom(p,q), (A (SAtom({sname="H"},_))) :: Number n :: l -> merge_atoms_inside id s (SAtom({p with shydrogens = int_of_string n},q),l)
  | SAtom(p,q), (A (SAtom({sname="H"},_))) :: l -> merge_atoms_inside id s (SAtom({p with shydrogens = 1},q),l)
  | SAtom(p,q), Symb "+" :: Number n :: l -> merge_atoms_inside id s (SAtom({p with scharge = int_of_string n},q),l)
  | SAtom(p,q), Symb "+" :: l -> merge_atoms_inside id s (SAtom({p with scharge = 1},q),l)
  | SAtom(p,q), Symb "-" :: Number n :: l -> merge_atoms_inside id s (SAtom({p with scharge = -(int_of_string n)},q),l)
  | SAtom(p,q), Symb "-" :: l -> merge_atoms_inside id s (SAtom({p with scharge = -1},q),l)
  | SAtom(p,q), Symb "@" :: l -> merge_atoms_inside id s (SAtom({p with schirality = AntiClockwise},q),l)
  | SAtom(p,q), Symb "@@" :: l -> merge_atoms_inside id s (SAtom({p with schirality = Clockwise},q),l)
  | SAtom(p,q), Symb ":" :: Number n :: RB :: l -> SAtom({p with sid = int_of_string n},q),l
  | SAtom(p,q), l -> failwith (s ^ " " ^ String.concat " :: " (Xlist.map l diagnostic_string_of_symbol))(*; SAtom(p,q),[]*)
  | _,_ -> failwith "merge_atoms_inside"

let rec merge_atoms id t = function
    A (SAtom(p,q)) :: l ->  A (SAtom(p,q)) :: (merge_atoms id t l)
  | Symb s :: l -> Symb s :: (merge_atoms id t l)
  | Number s :: l -> Number s :: (merge_atoms id t l)
  | LP :: l -> LP :: (merge_atoms id t l)
  | RP :: l -> RP :: (merge_atoms id t l)
  | LB :: A (SAtom(p,q)) :: l -> let a,l = merge_atoms_inside id t (SAtom(p,q),l) in (A a) :: (merge_atoms id t l)
  | [] -> []
  | l -> failwith (String.concat " :: " (Xlist.map l diagnostic_string_of_symbol))(*; []*)

let rec pointer_alpha_conversion map v = function
    Number n :: l ->
      (try
        let w = StringMap.find map n in
        Number (n ^ "." ^ w) :: (pointer_alpha_conversion (StringMap.remove map n) v l)
      with Not_found ->
        let w = string_of_int v in
        Number (n ^ "." ^ w) :: (pointer_alpha_conversion (StringMap.add map n w) (v+1) l))
  | s :: l -> s :: (pointer_alpha_conversion map v l)
  | [] -> []

let rec merge_pointers id t2 = function
    A (SAtom(p,q)) :: Number m :: l -> merge_pointers id t2 (A (SAtom(p,q @ [Default,SLink m])) :: l)
  | A (SAtom(p,q)) :: Symb "-" :: Number m :: l -> merge_pointers id t2 (A (SAtom(p,q @ [Single,SLink m])) :: l)
  | A (SAtom(p,q)) :: Symb "=" :: Number m :: l -> merge_pointers id t2 (A (SAtom(p,q @ [Double,SLink m])) :: l)
  | A (SAtom(p,q)) :: Symb "#" :: Number m :: l -> merge_pointers id t2 (A (SAtom(p,q @ [Triple,SLink m])) :: l)
  | A (SAtom(p,q)) :: Symb "/" :: Number m :: l -> merge_pointers id t2 (A (SAtom(p,q @ [Slash,SLink m])) :: l)
  | A (SAtom(p,q)) :: Symb "\\" :: Number m :: l -> merge_pointers id t2 (A (SAtom(p,q @ [Backslash,SLink m])) :: l)
  | A (SAtom(p,q)) :: l -> A (SAtom(p,q)) :: (merge_pointers id t2 l)
  | Symb s :: l -> Symb s :: (merge_pointers id t2 l)
  | LP :: l -> LP :: (merge_pointers id t2 l)
  | RP :: l -> RP :: (merge_pointers id t2 l)
  | [] -> []
  | l -> failwith (t2 ^ " " ^ String.concat " :: " (Xlist.map l diagnostic_string_of_symbol_simple))(*; []*)

let rec find_chain k rev = function
    A (SAtom(p,q)) :: l -> find_chain k (A (SAtom(p,q)) :: rev) l
  | Symb s :: l -> find_chain k (Symb s :: rev) l
  | LP :: l -> find_chain (k+1) (LP :: rev) l
  | RP :: l ->
      if k = 1 then List.rev rev, l else
      find_chain (k-1) (RP :: rev) l
  | l -> failwith ("find_chain" ^ String.concat " :: " (Xlist.map l diagnostic_string_of_symbol_simple))

let rec merge_connections id t2 = function
    [A (SAtom(p,q))] -> Default, SAtom(p,q)
  | A (SAtom(p,q)) :: LP :: l ->
      let x,y = find_chain 1 [] l in
      merge_connections id t2 (A (SAtom(p,q @ [merge_connections id t2 x])) :: y)
  | A (SAtom(p,q)) :: l -> Default, SAtom(p,q @ [merge_connections id t2 l])
  | Symb "-" :: A (SAtom(p,q)) :: l -> Single, snd (merge_connections id t2 (A (SAtom(p,q)) :: l))
  | Symb "=" :: A (SAtom(p,q)) :: l -> Double, snd (merge_connections id t2 (A (SAtom(p,q)) :: l))
  | Symb "#" :: A (SAtom(p,q)) :: l -> Triple, snd (merge_connections id t2 (A (SAtom(p,q)) :: l))
  | Symb "/" :: A (SAtom(p,q)) :: l -> Slash, snd (merge_connections id t2 (A (SAtom(p,q)) :: l))
  | Symb "\\" :: A (SAtom(p,q)) :: l -> Backslash, snd (merge_connections id t2 (A (SAtom(p,q)) :: l))
  | l -> failwith (t2 ^ " " ^ String.concat " :: " (Xlist.map l diagnostic_string_of_symbol_simple))(*; Default, SLink ""*)

let is_error = function
    ["0"; "."; a; b] -> true
  | _ -> false

let parse_smile_molecule id s =
  let l = Str.split (Str.regexp "") s in
  if is_error l then SLink("invalid smile " ^ s ^ " in record:") else
  let l = merge_atom_names id s l in
  let l = merge_numbers id s l in
  let l = merge_atoms id s l in
  let l = pointer_alpha_conversion StringMap.empty 0 l in
  let l = merge_pointers id s l in
  let t = snd (merge_connections id s l) in
  t

let parse_smile_reaction s =
  match Str.split (Str.regexp ">>") s with
    [a;b] -> Str.split (Str.regexp "\\.") a, Str.split (Str.regexp "\\.") b
  | _ -> failwith ("parse_smile_reaction: " ^ s)

(***********************************************************************************************)
(* translacja ze smilesu jako drzewa do smilesu jako napisu *)

let escape_string s =
  Int.fold 0 (String.length s - 1) "" (fun t i ->
    match String.sub s i 1 with
       "_" -> t ^ "\\_"
     | "#" -> t ^ "\\#"
     | "%" -> t ^ "\\%"
     | "\"" -> t ^ "''"
     | ">" -> t ^ "{}$>${}"
     | "\\" -> t ^ "{}$\\backslash${}"
     | c -> t ^ c)

let string_of_chirality = function
    AntiClockwise -> "@"
  | Clockwise -> "@@"
  | Unspecified -> ""

let string_of_bond = function
    Default -> ""
  | Single -> "-"
  | Slash -> "/"
  | Backslash -> "\\"
  | Double -> "="
  | Triple -> "#"
  | Aromatic -> ""

let latex_color = function
    0 -> "black"
  | 1 -> "red"
  | 2 -> "green"
  | 3 -> "blue"
  | 4 -> "magenta"
  | 5 -> "cyan"
  | 6 -> "Brown"
  | 7 -> "ForestGreen"
  | 8 -> "CarnationPink"
  | 9 -> "RoyalPurple"
  | 10 -> "LimeGreen"
  | 11 -> "OrangeRed"
  | 12 -> "Orange"
  | 13 -> "Sepia"
  | 14 -> "Rhodamine"
  | 15 -> "MidnightBlue"
  | 16 -> "Emerald"
  | 17 -> "BurntOrange"
  | 18 -> "Periwinkle"
  | 19 -> "Gray"
  | 20 -> "BrickRed"
  | 21 -> "OliveGreen"
  | 22 -> "SeaGreen"
  | 23 -> "YellowOrange"
  | 24 -> "Orchid"
  | 25 -> "BlueViolet"
  | 26 -> "PineGreen"
  | 27 -> "JungleGreen"
  | 28 -> "RawSienna"
  | 29 -> "Violet"
  | _ -> "yellow"

let gv_color = function
    0 -> "black"
  | 1 -> "red"
  | 2 -> "green"
  | 3 -> "blue"
  | 4 -> "magenta"
  | 5 -> "cyan"
  | 6 -> "Brown"
  | 7 -> "ForestGreen"
  | 8 -> "deeppink"
  | 9 -> "purple"
  | 10 -> "LimeGreen"
  | 11 -> "OrangeRed"
  | 12 -> "Orange"
  | 13 -> "thistle3"
  | 14 -> "violet"
  | 15 -> "MidnightBlue"
  | 16 -> "mediumslateblue"
  | 17 -> "orangered4"
  | 18 -> "slateblue3"
  | 19 -> "Gray"
  | 20 -> "indianred1"
  | 21 -> "khaki3"
  | 22 -> "SeaGreen"
  | 23 -> "deepskyblue2"
  | 24 -> "Orchid"
  | 25 -> "BlueViolet"
  | 26 -> "deeppink3"
  | 27 -> "mediumturquoise"
  | 28 -> "hotpink2"
  | 29 -> "maroon3"
  | _ -> "yellow"

type id_opts = NoId | StdId | StdLabel

let string_of_atom_bracets opts labels p =
  "[" ^
  (if p.svalence = Aromatic then String.lowercase_ascii  p.sname else p.sname) ^
  string_of_chirality p.schirality ^
  (match p.shydrogens with
    -1 -> ""
  | 0 -> ""
  | 1 -> "H"
  | n -> "H" ^ string_of_int n) ^
  (match p.scharge with
    0 -> ""
  | 1 -> "+"
  | -1 -> "-"
  | n -> if n > 0 then "+" ^ string_of_int n else "-" ^ string_of_int n) ^
  (if opts = StdId && p.sid <> 0 then ":" ^ string_of_int p.sid else "") ^
  (if opts = StdLabel && Labels.mem labels p.sid then ":" ^ string_of_int (Labels.get labels p.sid) else "") ^
  "]"

let string_of_atom opts labels p =
  if opts = StdLabel && Labels.mem labels p.sid then string_of_atom_bracets opts labels p else
  if opts = StdId && p.sid <> 0 then string_of_atom_bracets opts labels p else
  if p.scharge = 0 && p.schirality = Unspecified  && p.shydrogens = -1 then
    if p.svalence = Aromatic then String.lowercase_ascii  p.sname else p.sname
  else string_of_atom_bracets opts labels p

let rec string_of_tree opts labels = function
    SAtom(p,q) ->
      let s = string_of_atom opts labels p in
      s ^ string_of_tree_list opts labels q
  | SLink s ->
      let s = List.hd (Str.split (Str.regexp "\\.") s) in
      let s = if String.length s > 1 then "%" ^ s else s in
      s

and string_of_tree_list opts labels = function
    [] -> ""
  | [b,x] -> string_of_bond b ^ string_of_tree opts labels x
  | (b,SLink s) :: q -> string_of_bond b ^ string_of_tree opts labels (SLink s) ^ string_of_tree_list opts labels q
  | (b,x) :: q -> "(" ^ string_of_bond b ^ string_of_tree opts labels x ^ ")" ^ string_of_tree_list opts labels q

let string_of_tree_no_id t = string_of_tree NoId (Labels.make 1) t

let string_of_tree_std_id t = string_of_tree StdId (Labels.make 1) t

(***********************************************************************************************)
(* wykrywanie cykli *)

let rec get_links links = function
  | Atom(_,l) -> Xlist.fold l links (fun links (_,a) -> get_links links a)
  | Link t -> StringSet.add links t

let rec get_link_paths path map = function
  | Atom(p,l) -> Xlist.fold l map (fun map (_,a) -> get_link_paths (p.id :: path) map a)
  | Link t -> StringMap.add_inc map t [List.rev path] (fun l -> (List.rev path) :: l)

let rec extract_cycles_rec t s = function
    n1 :: m1 :: p1, n2 :: m2 :: p2 ->
      if n1 <> n2 then failwith "extract_cycles_rec 1a" else
      if m1 = m2 then extract_cycles_rec t s (m1 :: p1,m2 :: p2) else
      [n1;m1;m2] @ p1 @ p2
  | [n1], n2 :: p2 ->
      if n1 <> n2 then failwith "extract_cycles_rec 1b" else
      [n1] @ p2
  | n1 :: p1, [n2] ->
      if n1 <> n2 then failwith "extract_cycles_rec 1c" else
      [n2] @ p1
  | p1,p2 ->
(*      Printf.printf "extract_cycles: bad paths [%s] [%s] to link %s in %s\n%!"
        (String.concat ";" (Xlist.map p1 string_of_int)) (String.concat ";" (Xlist.map p2 string_of_int)) t s; *)
      failwith "extract_cycles_rec 2"

let extract_cycles s paths =
  StringMap.fold paths [] (fun cycles t -> function
      [p1;p2] -> (extract_cycles_rec t s (p1,p2)) :: cycles
    | _ -> (*Printf.printf "extract_cycles: bad link %s in %s\n%!" t s;*) cycles)

let extract_labelled_cycles s paths =
  StringMap.mapi paths (fun t -> function
      [p1;p2] -> (extract_cycles_rec t s (p1,p2))
    | _ -> (*Printf.printf "extract_labelled_cycles: bad link %s in %s\n%!" t s;*) [])

let transitive_closure map =
  let f = ref true in
  let r = ref map in
  while !f do
    f := false;
    r := StringMap.fold (!r) (!r) (fun map k set ->
      let set2 = StringSet.fold set set (fun set2 v ->
        StringSet.union set2 (StringMap.find map v)) in
      if StringSet.size set2 > StringSet.size set then f := true;
      StringMap.add map k set2)
  done;
  !r

let cycle_union cycles =
  let map = StringMap.fold cycles IntMap.empty (fun map t l ->
    Xlist.fold l map (fun map n ->
      IntMap.add_inc map n [t] (fun l -> t :: l))) in
  let labels = IntMap.fold map StringMap.empty (fun labels _ l ->
    let set = Xlist.fold l StringSet.empty StringSet.add in
    Xlist.fold l labels (fun labels v ->
      StringMap.add_inc labels v set (fun set2 -> StringSet.union set set2))) in
  let labels = transitive_closure labels in
  let labels,_ = StringMap.fold labels ([],StringSet.empty) (fun (l,sum) v set ->
    if StringSet.mem sum v then l,sum else
    set :: l, StringSet.union sum set) in
  Xlist.map labels (fun set ->
    set, StringSet.fold set IntSet.empty (fun set v ->
      Xlist.fold (StringMap.find cycles v) set IntSet.add))

let find_cycles t =
    let paths = get_link_paths [] StringMap.empty t in
    let cycles = extract_labelled_cycles "" paths in
    let cycles = cycle_union cycles in
    cycles

(***********************************************************************************************)
(* nadawanie id oraz translacja ze smilesu jako drzewa do cząsteczki jako drzewa *)

let int_of_bond1 = function
    Default -> 1
  | Single -> 1
  | Slash -> 1
  | Backslash -> 1
  | Double -> 2
  | Triple -> 3
  | Aromatic -> failwith "int_of_bond"

let int_of_bond = function
    Default -> 2
  | Single -> 2
  | Slash -> 2
  | Backslash -> 2
  | Double -> 4
  | Triple -> 6
  | Aromatic -> 3

let lowest_normal_valence n = function
    "B" -> if n > 3 then failwith "lowest_normal_valence B" else 3
  | "C" -> if n > 4 then failwith ("lowest_normal_valence C " ^ string_of_int n) else 4
  | "N" -> if n > 5 then failwith "lowest_normal_valence N" else if n > 3 then 5 else 3
  | "O" -> if n > 2 then failwith "lowest_normal_valence O" else 2
  | "P" -> if n > 5 then failwith "lowest_normal_valence P" else if n > 3 then 5 else 3
  | "S" -> if n > 6 then failwith "lowest_normal_valence S" else if n > 4 then 6 else if n > 2 then 4 else 2
  | "F" -> if n > 1 then failwith "lowest_normal_valence F" else 1
  | "Cl" -> if n > 1 then failwith "lowest_normal_valence Cl" else 1
  | "Br" -> if n > 1 then failwith "lowest_normal_valence Br" else 1
  | "I" -> if n > 1 then n(*failwith "lowest_normal_valence I"*) else 1
  | _ -> failwith "lowest_normal_valence"

let rec calculate_no_hydrogens n0 = function
    SAtom(p,q) ->
      let no_hydrogens =
        if p.shydrogens = -1 then
          match p.svalence with
            Alifatic ->
              let n = Xlist.fold q n0 (fun n (b,_) -> int_of_bond1 b + n) in
              let v = lowest_normal_valence n p.sname in
              v-n
          | Aromatic ->
              if p.sname = "C" then if Xlist.size q = 3 then 0 else 1
              else 0
        else p.shydrogens in
      let q = Xlist.map q (fun (b,a) -> b, calculate_no_hydrogens (int_of_bond1 b) a) in
      SAtom({p with shydrogens=no_hydrogens},q)
  | SLink t -> SLink t

let rec manage_bonds links = function
    SAtom(p,q) ->
      let q = Xlist.map q (fun (b,a) ->
        if b <> Default then b,a else
        if p.svalence <> Aromatic then Single,a else
        match a with
          SAtom({svalence=Aromatic},_) -> (*Default*)Aromatic,a
        | SAtom({svalence=Alifatic},_) -> Single,a
        | SLink t -> if StringSet.mem links t then (*Default*)Aromatic,a else Single,a) in
      SAtom(p,Xlist.map q (fun (b,a) -> b, manage_bonds links a))
  | SLink t -> SLink t

let rec find_aromatic_links aromatic links = function
    SAtom(p,q) ->
      Xlist.fold q links (fun links (b,a) ->
        let aromatic = p.svalence = Aromatic && b = Default in
        find_aromatic_links aromatic links a)
  | SLink t -> if aromatic then StringQMap.add links t else links

(*let aromatic_schemata = Xlist.fold [
  "C",[[1;1;2];[1;2;1];[2;1;1];[1;2];[2;1]];
  "N",[[1;2];[2;1]];
  "NH",[[1;1]];
  "S",[[1;1]];
  "O",[[1;1]];
  ] StringMap.empty (fun map (k,v) -> StringMap.add map k v)

let rec manage_aromaticity_rec links input = function
    SAtom(p,q) ->
      let schemata =
        if p.svalence = Aromatic then
          if p.sname = "N" && p.shydrogens = 1 then StringMap.find aromatic_schemata "NH" else
          try StringMap.find aromatic_schemata p.sname with Not_found -> failwith ("manage_aromaticity_rec: aromatic " ^ p.sname)
        else [[]] in
      let schemata = if input = 0 then schemata else
        Xlist.fold schemata [] (fun schemata -> function
          x :: schema -> if x = input then schema :: schemata else schemata
        | [] -> schemata) in
      Xlist.fold schemata [] (fun found schema ->
        let l = manage_aromaticity_list schema links [] q in
        Xlist.fold l found (fun found (q,links) ->
          (SAtom({p with svalence = Alifatic},q), links) :: found))
  | SLink t ->
      if StringMap.mem links t then
        if StringMap.find links t = input then [SLink t, links] else []
      else [SLink t, StringMap.add links t input]

and manage_aromaticity_list schema links q_rev = function
    [] -> if schema = [] then [List.rev q_rev,links] else []
  | (b,a) :: q ->
          if b = Default && schema = [] then [] else
          let b,input,schema = if b <> Default then b,0,schema else
            if List.hd schema = 1 then Single,1,List.tl schema else Double,2,List.tl schema in
          let l = manage_aromaticity_rec links input a in
          Xlist.fold l [] (fun l (a,links) ->
            (manage_aromaticity_list schema links ((b,a) :: q_rev) q) @ l)*)

let manage_aromaticity tree =
  let links = find_aromatic_links false StringQMap.empty tree in
  let links = StringQMap.fold links StringSet.empty (fun links t v ->
    if v = 2 then StringSet.add links t else links) in
  let tree = manage_bonds links tree in
  tree
(*  match manage_aromaticity_rec StringMap.empty 0 tree with
    [] -> failwith "manage_aromaticity"
(*   | (tree,_) :: _ -> tree *)
  | l -> Xlist.rev_map l fst*)

let rec sort_links reva revl = function
    [] -> List.rev revl @ List.rev reva
  | (b,Atom(p,q)) :: l -> sort_links ((b,Atom(p,q)) :: reva) revl l
  | (b,Link t) :: l -> sort_links reva ((b,Link t) :: revl) l

(* FIXME: problem z O podłączonym podwójnym wiązaniem z C aromatycznym *)
let rec atom_tree_of_smile_tree = function
    SAtom(p,q) ->
      let hydrogens = Int.fold 1 p.shydrogens [] (fun hydrogens _ ->
            (Single,Atom({name="H"; charge=0; chirality=Unspecified; id=0; hydrogens=[]; fluors=[]},[])) :: hydrogens) in
      let q = Xlist.map q (fun (b,a) -> b, atom_tree_of_smile_tree a) in
      Atom({name=p.sname; charge=p.scharge; chirality=p.schirality; id=p.sid; hydrogens=[]; fluors=[]},(*sort_links [] []*) (hydrogens @ q)) (* FIXME: sort_links zaburza chiralność *)
  | SLink t -> Link t

let make_molecule stoi s =
  let a = atom_tree_of_smile_tree (manage_aromaticity (calculate_no_hydrogens 0 (parse_smile_molecule "" s))) in
  {smiles=s; smile_tree=SLink "ni"; tree=a; ids=IntSet.empty;
    hf_ids=IntSet.empty; stoi=stoi; obligatory=false; cycles=[]}

(**
let rec assign_unique_ids n = function
    SAtom(p,q) ->
     let p = {p with sid=n} in
     let q,n = Xlist.fold q ([],n+1) (fun (q,n) (b,a) ->
       let a,n = assign_unique_ids n a in (b,a) :: q, n) in
     SAtom(p,List.rev q),n
  | SLink t -> SLink t,n

**)

let rec assign_unique_ids n = function
  | Atom(p,q) ->
     let p = {p with id=n} in
     let q,n = Xlist.fold q ([],n+1) (fun (q,n) (b,a) ->
       let a,n = assign_unique_ids n a in (b,a) :: q, n) in
     Atom(p,List.rev q),n
  | Link t -> Link t,n

let rec get_tree_ids set = function
    Atom(p,q) ->
      let set = IntSet.add set p.id in
      Xlist.fold q set (fun set (_,a) -> get_tree_ids set a)
  | Link _ -> set

let assign_unique_ids_list n l =
  let l,n = Xlist.fold l ([],n) (fun (l,n) m ->
    let a,n = assign_unique_ids n m.tree in {m with tree=a; ids=get_tree_ids IntSet.empty a} :: l, n) in
  List.rev l, n

let make_molecule_pattern s =
  let a = atom_tree_of_smile_tree (manage_aromaticity (calculate_no_hydrogens 0 (parse_smile_molecule "" s))) in
  let ids = get_tree_ids IntSet.empty a in
  {smiles=s; smile_tree=SLink "ni"; tree=a; ids=ids;
    hf_ids=IntSet.empty; stoi=false; obligatory=false; cycles=[]}

let bond_of_int = function
    2 -> Single
  | 4 -> Double
  | 6 -> Triple
  | 3 -> Aromatic
  | n -> failwith ("bond_of_int: " ^ string_of_int n)

(**
let rec atom_tree_of_smile_tree_simple = function
    SAtom(p,q) ->
      let q = Xlist.map q (fun (b,a) -> b, atom_tree_of_smile_tree_simple a) in
      Atom({name=p.sname; charge=p.scharge; id=p.sid; hydrogens=[]; fluors=[]},q)
  | SLink t -> Link t

let rec add_unique_ids n = function
    Atom(p,q) ->
     let p,n = if p.id = 0 then {p with id=n}, n+1 else p,n in
     let q,n = Xlist.fold q ([],n) (fun (q,n) (b,a) ->
       let a,n = add_unique_ids n a in (b,a) :: q, n) in
     Atom(p,List.rev q),n
  | Link t -> Link t,n

let assign_unique_ids_and_translate_to_tree n t =
  let t,n = assign_unique_ids n t in
  let t2 = atom_tree_of_smile_tree 0 t in
  let t2,n = add_unique_ids n t2 in
  t,t2,n

let assign_unique_ids_and_translate_to_tree_list n l =
  let l,l2,n = Xlist.fold l ([],[],n) (fun (l,l2,n) a ->
    let t,t2,n = assign_unique_ids_and_translate_to_tree n a in t :: l, t2 :: l2, n) in
  List.rev l, List.rev l2, n


(*let make_molecule stoi (rev,n) s =
  let al = Xlist.rev_map (manage_aromaticity (parse_smile_molecule "" s)) (atom_tree_of_smile_tree 0)  in
  let anl = Xlist.rev_map al (assign_unique_ids_vis n) in
  let n = snd (List.hd anl) in
  (Xlist.rev_map anl (fun (a,_) -> {smiles=s; smile_tree=SLink "ni"; tree=a; ids=get_tree_ids IntSet.empty a;
    hf_ids=IntSet.empty; stoi=stoi; cycles=find_cycles a})) :: rev, n*)
let make_molecule stoi (rev,n) s =
  let a = atom_tree_of_smile_tree 0 (manage_aromaticity (calculate_no_hydrogens 0 (parse_smile_molecule "" s)))  in
  let a,n = assign_unique_ids_vis n a in
  {smiles=s; smile_tree=SLink "ni"; tree=a; ids=get_tree_ids IntSet.empty a;
    hf_ids=IntSet.empty; stoi=stoi; cycles=find_cycles a} :: rev, n
    **)
(***********************************************************************************************)
(* translacja do xml *)

let rec smile_tree_of_atom_tree = function
    Atom(p,q) ->
      let arom = Xlist.fold q false (fun arom (b,_) -> if b = Aromatic then true else arom) in
      let q = Xlist.map q (fun (b,a) -> (if b = Aromatic then Default else b), smile_tree_of_atom_tree a) in
      SAtom({sname=p.name; scharge=p.charge; schirality=p.chirality; sid=p.id; svalence=if arom then Aromatic else Alifatic; shydrogens=0},q)
  | Link t -> SLink t

let solution_to_xml r labels center_bonds reid =
        Xml.Element("solution",[],[
          Xml.Element("atoms",[],List.rev (Int.fold 1 (r.reaction_size-1) [] (fun l i ->
            Xml.Element("atom",["id",string_of_int i;"visible",if Labels.get labels i = -1 then "f" else "t"],
              (if Labels.get labels i <> -1 then [Xml.Element("colour",[],[Xml.PCData (latex_color (Labels.get labels i))])] else []) @
              (let id = try IntMap.find reid i with Not_found -> i in
               if Labels.get labels i <> -1 && id <> 0 then [Xml.Element("label",[],[Xml.PCData (string_of_int id)])] else [])) :: l)));
            Xml.Element("bonds",[],
              List.rev (IntMap.fold center_bonds [] (fun l i map ->
                IntMap.fold map l (fun l j label ->
                  Xml.Element("bond",["id1",string_of_int i;"id2",string_of_int j],[
                    Xml.Element("colour",[],[Xml.PCData (latex_color label)])]) :: l))))])

let solution_to_xml2 (r,labels,_,center_bonds,reid) =
      let xml_reactants = String.concat "." (Xlist.map r.reactants (fun m ->
        string_of_tree_std_id (smile_tree_of_atom_tree m.tree))) in
      let xml_products = String.concat "." (Xlist.map r.products (fun m ->
        string_of_tree_std_id (smile_tree_of_atom_tree m.tree))) in
      Xml.Element("solution",[],[
          Xml.Element("reactants",[],[Xml.PCData xml_reactants]);
          Xml.Element("products",[],[Xml.PCData xml_products]);
          Xml.Element("atoms",[],List.rev (Int.fold 1 (r.reaction_size-1) [] (fun l i ->
            Xml.Element("atom",["id",string_of_int i;"visible",if Labels.get labels i = -1 then "f" else "t"],
              (if Labels.get labels i <> -1 then [Xml.Element("colour",[],[Xml.PCData (latex_color (Labels.get labels i))])] else []) @
              (let id = try IntMap.find reid i with Not_found -> i in
               if Labels.get labels i <> -1 && id <> 0 then [Xml.Element("label",[],[Xml.PCData (string_of_int id)])] else [])) :: l)));
            Xml.Element("bonds",[],
              List.rev (IntMap.fold center_bonds [] (fun l i map ->
                IntMap.fold map l (fun l j label ->
                  Xml.Element("bond",["id1",string_of_int i;"id2",string_of_int j],[
                    Xml.Element("colour",[],[Xml.PCData (latex_color label)])]) :: l))))])

let reaction_to_xml r msg labels center_bonds reid =
      let xml_reactants = String.concat "." (Xlist.map r.reactants (fun m ->
        string_of_tree_std_id (smile_tree_of_atom_tree m.tree))) in
      let xml_products = String.concat "." (Xlist.map r.products (fun m ->
        string_of_tree_std_id (smile_tree_of_atom_tree m.tree))) in
      let xml_solution = solution_to_xml r labels center_bonds reid in
      let xml =
            Xml.Element("reaction",[],[
              Xml.Element("reactants",[],[Xml.PCData xml_reactants]);
              Xml.Element("products",[],[Xml.PCData xml_products]);
              Xml.Element("messages",[],Xlist.map msg (fun s -> Xml.Element("message",[],[Xml.PCData s])));
              Xml.Element("solutions",[],[xml_solution])]) in
      xml

let solutions_to_xml reaction_smile messages solutions =
  let solutions = Xlist.map solutions solution_to_xml2 in
  let xml =
    Xml.Element("reaction",[],[
      Xml.Element("smiles",[],[Xml.PCData reaction_smile]);
      Xml.Element("messages",[],Xlist.map messages (fun s -> Xml.Element("message",[],[Xml.PCData s])));
      Xml.Element("solutions",[],solutions)]) in
  xml

(* FIXME: niepoprawnie wyświetla fluory po ich zwinięciu *)
let reaction_to_gv file min_i max_i r msg labels center_bonds reid =
  Printf.fprintf file "graph G {\n  node [shape=plaintext,margin=0,fontsize=30]\n";
  Int.iter min_i (max_i-1) (fun i ->
    if Labels.get labels i = -1 then () else (
      let p,q = r.graph.(i) in
      let id = try IntMap.find reid i with Not_found -> i in
      if id = 0 then Printf.fprintf file "  %d [label=<%s>,fontcolor=%s];\n" p.id p.name (gv_color (Labels.get labels p.id))
      else Printf.fprintf file "  %d [label=<<SUP>%d</SUP>%s>,fontcolor=%s];\n" p.id id p.name (gv_color (Labels.get labels p.id));
      Xlist.iter p.fluors (fun j -> Printf.fprintf file "  %d -- %d;\n" p.id j);
      Xlist.iter q (fun (b,a) ->
        if p.id < a.id && Labels.get labels a.id <> -1 then
          let color = gv_color (try IntMap.find (IntMap.find center_bonds p.id) a.id with Not_found -> 0) in
          match b with
          | Double -> Printf.fprintf file "  %d -- %d [color=%s];\n" p.id a.id color; Printf.fprintf file "  %d -- %d [color=%s];\n" p.id a.id color
          | Triple -> Printf.fprintf file "  %d -- %d [color=%s];\n" p.id a.id color; Printf.fprintf file "  %d -- %d [color=%s];\n" p.id a.id color; Printf.fprintf file "  %d -- %d [color=%s];\n" p.id a.id color
          | Aromatic -> Printf.fprintf file "  %d -- %d [color=%s];\n" p.id a.id color; Printf.fprintf file "  %d -- %d [style=dashed,color=%s];\n" p.id a.id color
          | _ -> Printf.fprintf file "  %d -- %d [color=%s];\n" p.id a.id color)));
(*       Printf.fprintf file "  %d -> %d  [%s]\n" p.id a.id "x")); *)
  Printf.fprintf file "  }\n"

(***********************************************************************************************)
(* translacja cząsteczki jako drzewa do cząsteczki jako grafu *)

let rec make_links_map links p = function
    _,Atom(p,l) -> Xlist.fold l links (fun links a -> make_links_map links p a)
  | b,Link t -> StringMap.add_inc links t (b,p,[]) (function (b2,p2,[]) -> b2,p2,[b,p] | _ -> failwith "make_links_map")

let add_edge graph from_node to_node bound =
  if (fst (graph.(from_node.id))).id = 0 then
    graph.(from_node.id) <- from_node,[bound,to_node]
  else
    let _,l = graph.(from_node.id) in
    graph.(from_node.id) <- from_node,(bound,to_node) :: l

let rec make_atom_graph_rec graph q = function
    b,Atom(p,l) ->
      add_edge graph p q b;
      add_edge graph q p b;
      Xlist.iter l (fun a -> make_atom_graph_rec graph p a)
  | _,Link t -> ()

let match_bounds t = function
    Default,Single ->  Single
  | Single,Default -> Single
  | Single,Slash -> Single (* FIXME: tu powinien byc blad *)
  | Single,Backslash -> Single
  | Slash,Single -> Single
  | Backslash,Single -> Single
  | Default,x -> x (* FIXME: tu powinien byc blad *)
  | x,Default -> x
  | Double,Single -> Double
  | Single,Double -> Double
  | Triple,Single -> Double
  | Single,Triple -> Double
  | _ -> failwith ("Inconsistent bonds near link " ^ (List.hd (Xstring.split "\\." t)))

let make_atom_graph size l =
  let graph = Array.make size ({name=""; charge=0; id=0; hydrogens=[]; chirality=Unspecified; fluors=[]},[]) in
  Xlist.iter l (function
      Atom(p,l) ->
          (* print_endline (string_of_tree_no_id (smile_tree_of_atom_tree (Atom(p,l)))); *)
          graph.(p.id) <- p,[];
          Xlist.iter l (fun a -> make_atom_graph_rec graph p a);
          let links = Xlist.fold l StringMap.empty (fun links a -> make_links_map links p a) in
          StringMap.iter links (fun t -> function
              (b,p,[b2,q]) ->
                 let b = if b = b2 then b else match_bounds t (b,b2) in
                 add_edge graph p q b;
                 add_edge graph q p b
            | _ -> failwith "make_atom_graph 3a")
    | Link t -> ()(*failwith ("make_atom_graph 1: " ^ t)*));
  graph

(***********************************************************************************************)
(* wizualizacja w latexu *)

let a0poster_header papersize =
  String.concat "\n" [
    "\\documentclass{article}";
(*     "\\documentclass[portrait,final]{a0poster}"; *)
    "\\usepackage[" ^ papersize ^ "paper, left=2.5cm, right=2.5cm, top=3.5cm, bottom=3.5cm, headsep=1.2cm]{geometry}";
(*     "\\newgeometry{tmargin=3cm, bmargin=3cm, lmargin=0.5cm, rmargin=0.5cm}"; *)
(* "\\usepackage{a4wide}"; *)
(*     "\\usepackage{color}"; *)
    "\\usepackage[usenames,dvipsnames]{xcolor}";
    "\\usepackage{amsmath}";
    "\\usepackage{amssymb}";
(*     "\\usepackage[T1]{fontenc}"; *)
(*     "\\usepackage[utf8]{inputenc}"; *)
(*     "\\usepackage[polish]{babel}"; *)
(*     "\\usepackage{tikz}";  *)
(*     "\\usetikzlibrary{conceptgraph}";  *)
    "\\usepackage{chemfig}";
    "\\renewcommand*\\printatom[1]{\\ensuremath{\\mathsf{#1}}}";
    "\\parindent 0pt";
    "\\parskip 4pt";
    "\\begin{document}\n\n"]

let trailer = "\\end{document}"
(*
let latex_chemfig_of_bond = function
    Default -> "-"
  | Single -> "-"
  | Slash -> "-"
  | Backslash -> "-"
  | Double -> "="
  | Triple -> "~"

let latex_chemfig_of_charge = function
    0 -> ""
  | 1 -> "^+"
  | -1 -> "^-"
  | n -> if n > 0 then "^{" ^ string_of_int n ^ "+}" else "^{" ^ string_of_int (-n) ^ "-}"

(*let rec latex_chemfig_of_node labels angles = function
    b,Atom(p,l) -> Printf.sprintf "%s[:%d]{\\textcolor{%s}{%s_{%d}%s}}%s"
                    (latex_chemfig_of_bond b)
                    (try angles.(p.id) with _ -> Printf.printf "latex_chemfig_of_node %s:%d \n%!" p.name p.id; 0)
                    (latex_color (Labels.get labels p.id))
                    p.name p.id
                    (latex_chemfig_of_charge p.charge)
                    (latex_chemfig_of_list labels angles l)
  | b,Link s -> Printf.sprintf "?[%s,%d]" s (int_of_bond b)

and latex_chemfig_of_list labels angles = function
    [] -> ""
  | [x] -> latex_chemfig_of_node labels angles x
  | x :: l -> "(" ^ latex_chemfig_of_node labels angles x ^ ")" ^ latex_chemfig_of_list labels angles l

let latex_chemfig_of labels angles tree =
  "\\chemfig{" ^
  (match tree with
    Atom(p,l) -> Printf.sprintf "{\\textcolor{%s}{%s_{%d}%s}}%s"
                    (latex_color (Labels.get labels p.id))
                    p.name p.id
                    (latex_chemfig_of_charge p.charge)
                    (latex_chemfig_of_list labels angles l)
  | Link _ -> failwith "latex_chemfig_of")
  ^ "}"*)

let latex_chemfig_of_id reid id =
  let id = try IntMap.find reid id with Not_found -> id in
  if id = 0 then "" else "_{" ^ string_of_int id ^ "}"

let rec latex_chemfig_of_node id labels bond_labels link_bond_labels angles reid = function
    b,Atom(p,l) ->
       if Labels.get labels p.id = -1 then "" else
       Printf.sprintf "%s[:%d%s]{\\textcolor{%s}{%s%s%s}}%s"
                    (latex_chemfig_of_bond b)
                    angles.(p.id)
                    (try ",,,," ^ latex_color (IntMap.find (IntMap.find bond_labels (min p.id id)) (max p.id id)) with Not_found -> "")
                    (latex_color (Labels.get labels p.id))
                    p.name (latex_chemfig_of_id reid p.id) (latex_chemfig_of_charge p.charge)
                    (latex_chemfig_of_list p.id labels bond_labels link_bond_labels angles reid l)
  | b,Link s -> Printf.sprintf "?[%s,%d%s]" s (int_of_bond b) (try "," ^ latex_color (StringMap.find link_bond_labels s) with Not_found -> "")

and latex_chemfig_of_list id labels bond_labels link_bond_labels angles reid = function
    [] -> ""
  | [x] -> latex_chemfig_of_node id labels bond_labels link_bond_labels angles reid x
  | x :: l -> "(" ^ latex_chemfig_of_node id labels bond_labels link_bond_labels angles reid x ^ ")" ^ latex_chemfig_of_list id labels bond_labels link_bond_labels angles reid l

let rec latex_chemfig_of labels bond_labels angles reid = function
    Atom(p,l) ->
       if Labels.get labels p.id = -1 && Xlist.size l > 0 then (* FIXME: prawdopodobna przyczyna żółtych wodorów *)
(*          if Xlist.size l = 0 then "HH" else *)
         if Xlist.size l = 1 then latex_chemfig_of labels bond_labels angles reid (snd (List.hd l))
         else failwith "latex_chemfig_of 2" else
       let links = Xlist.fold l StringMap.empty (fun links a -> make_links_map links p a) in
       let link_bond_labels = StringMap.fold links StringMap.empty (fun link_bond_labels t -> function
            (_,p,[_,q]) -> (try StringMap.add link_bond_labels t (IntMap.find (IntMap.find bond_labels (min p.id q.id)) (max p.id q.id)) with Not_found -> link_bond_labels)
          | _ -> failwith "latex_chemfig_of 3") in
       Printf.sprintf "\\chemfig{{\\textcolor{%s}{%s%s%s}}%s}"
                    (latex_color (Labels.get labels p.id))
                    p.name (latex_chemfig_of_id reid p.id) (latex_chemfig_of_charge p.charge)
                    (latex_chemfig_of_list p.id labels bond_labels link_bond_labels angles reid l)
  | Link _ -> failwith "latex_chemfig_of 1"
*)
(***********************************************************************************************)

let rec calculate_stoi2 map = function
    Atom(p,l) -> Xlist.fold l (StringQMap.add map p.name) (fun map (_,a) -> calculate_stoi2 map a)
  | Link _ -> map

let calculate_stoi_list l =
  Xlist.fold l StringQMap.empty (fun qmap m -> calculate_stoi2 qmap m.tree)

let calculate_stoi reactants products =
  let reactant_stoi = calculate_stoi_list reactants in
  let product_stoi = calculate_stoi_list products in
  let stoi = StringQMap.fold reactant_stoi StringMap.empty (fun stoi s v -> StringMap.add stoi s (v,0)) in
  let stoi = StringQMap.fold product_stoi stoi (fun stoi s v -> StringMap.add_inc stoi s (0,v) (fun (w,_) -> w,v)) in
  let stoi = StringMap.fold stoi StringMap.empty (fun stoi name (v,w) -> if v = w then stoi else StringMap.add stoi name (v,w)) in
  stoi

let collapse_hydrogens_and_fluors graph = (* FIXME: sprawdzenie czy wiązanie jest pojedyncze i czy atom jest jednowartościowy *)
  Array.map (fun (p,l) ->
    if p.name = "H" || p.name = "F" then p,l else
    let p,l = Xlist.fold l (p,[]) (fun (p,l) (b,q) ->
      if q.name = "H" then {p with hydrogens=q.id :: p.hydrogens},l else
      if q.name = "F" then {p with fluors=q.id :: p.fluors},l else p,(b,q) :: l) in
    p,List.rev l) graph

let select_non_hf_atoms graph set =
  IntSet.fold set (IntSet.empty,IntSet.empty) (fun (ids,hf_ids) id ->
    if (fst graph.(id)).name = "H" || (fst graph.(id)).name = "F" then ids,IntSet.add hf_ids id else IntSet.add ids id,hf_ids)

let find_duplicates stoi reactants products messages =
  let core_stoi = List.sort compare (StringMap.fold stoi [] (fun l name (v,w) -> if name = "H" || name = "O" then l else (name,w - v) :: l)) in
  let duplicates = if core_stoi = [] || Xlist.fold core_stoi false (fun b (_,v) -> if v < 0 then true else b) then [] else
    Xlist.fold reactants [] (fun duplicates t ->
      let reactant_stoi = calculate_stoi2 StringQMap.empty t.tree in
      let reactant_core_stoi = List.sort compare (StringQMap.fold reactant_stoi [] (fun l name v -> if name = "H" || name = "O" then l else (name,v) :: l)) in
      if core_stoi = reactant_core_stoi then t :: duplicates else duplicates) in
  let messages = if Xlist.size duplicates > 1 then messages @ [Printf.sprintf "MULTIPLE CANDIDATES FOR DUPLICATES: %s" (escape_string (String.concat "." (Xlist.map duplicates (fun m -> m.smiles))))] else messages in
  (if duplicates = [] then [] else [List.hd duplicates]), messages

let add_waters stoi =
  if StringMap.size stoi <> 2 || not (StringMap.mem stoi "H") || not (StringMap.mem stoi "O") then [],[] else
  let hydrogens = let v,w = StringMap.find stoi "H" in v - w in
  let oxygens = let v,w = StringMap.find stoi "O" in v - w in
  if hydrogens = oxygens * 2 then
    if oxygens > 0 then [], Int.fold 1 oxygens [] (fun l _ -> make_molecule false "O" :: l)
    else Int.fold 1 (-oxygens) [] (fun l _ -> make_molecule false "O" :: l),[] else [],[]

let repair_stoi_reaction reactants products messages =
  let stoi = calculate_stoi reactants products in
  if StringMap.is_empty stoi then reactants,products,messages else
  let messages = if StringMap.is_empty stoi then messages else
    let l = List.sort compare (StringMap.fold stoi [] (fun l name (v,w) -> (Printf.sprintf "%s: %d$>>$%d" name v w) :: l)) in
    messages @ [Printf.sprintf "INVALID STOICHIOMETRY: %s" (escape_string (String.concat " " l))] in
  let added_reactants,messages = find_duplicates stoi reactants products messages in
  let stoi = calculate_stoi (added_reactants @ reactants) products in
  let added_reactant_waters, added_product_waters = add_waters stoi in
  let added_reactants = added_reactants @ added_reactant_waters in
  let added_products = added_product_waters in
  let messages = if added_reactants <> [] then messages @ [Printf.sprintf "added reactants: %s" (escape_string (String.concat "." (Xlist.map added_reactants (fun t -> t.smiles))))] else messages in
  let messages = if added_products <> [] then messages @ [Printf.sprintf "added products: %s" (escape_string (String.concat "." (Xlist.map added_products (fun t -> t.smiles))))] else messages in
  let reactants = reactants @ added_reactants in
  let products = products @ added_products in
  let stoi = calculate_stoi reactants products in
  let reactants,products = StringMap.fold stoi (reactants,products) (fun (reactants,products) name (v,w) ->
    let t = make_molecule true ("[" ^ name ^ "]") in
    if v < w then Int.fold 1 (w-v) reactants (fun reactants _ -> t :: reactants), products
    else reactants, Int.fold 1 (v-w) products (fun products _ -> t :: products)) in
  List.rev reactants,List.rev products,messages

let prepare_record_for_matching r messages =
(*     let r = find_aux_solvents r in *)
(*   printf "prepare_record_for_matching 1\n%!"; *)
  let reactant_smiles,product_smiles = parse_smile_reaction r.reaction_smile in
  let reactants,products = Pair.list_map (reactant_smiles,product_smiles) (make_molecule false) in
  let reactants,products,messages = repair_stoi_reaction reactants products messages in
  let reactants,ids_reactants = assign_unique_ids_list 1 reactants in
  let products,ids_products = assign_unique_ids_list ids_reactants products in
  let graph = make_atom_graph ids_products (Xlist.map (reactants @ products) (fun m -> m.tree)) in
  let graph = collapse_hydrogens_and_fluors graph in
  let reactants = Xlist.map reactants (fun m -> let ids,hf_ids = select_non_hf_atoms graph m.ids in {m with ids=ids; hf_ids=hf_ids}) in
  let products = Xlist.map products (fun m -> let ids,hf_ids = select_non_hf_atoms graph m.ids in {m with ids=ids; hf_ids=hf_ids}) in
  let reactant_ids = Xlist.fold reactants IntSet.empty (fun set m -> IntSet.union set m.ids) in
  let product_ids = Xlist.fold products IntSet.empty (fun set m -> IntSet.union set m.ids) in
(*       printf "prepare_record_for_matching 4b\n%!"; *)
(*    let (reactant_ids,one_atom_reactants),(product_ids,one_atom_products) = Pair.list_fold (reactant_trees,product_trees) ((reactant_ids,IntSet.empty),(product_ids,IntSet.empty)) (fun (setm,seto) t ->
      let id = is_one_atom_tree t in
      if id = -1 then setm, seto else IntSet.remove setm id, IntSet.add seto id) in*)
(*  let cycles = Xlist.fold (reactants @ products) IntSet.empty (fun cycles m ->
    Xlist.fold m.cycles cycles (fun cycles (_,set) -> IntSet.union set cycles)) in*)
(*       printf "prepare_record_for_matching 4c\n%!"; *)
  let empty_labels = Labels.make ids_products in
(*       printf "prepare_record_for_matching 4d\n%!"; *)
(*       let reactant_groups = CommonSubstructure.find_isomorphic_components graph (Collection.of_list (Xlist.map reactants (fun reactant -> Collection.of_list (IntSet.to_list reactant.ids)))) empty_labels in *)
(*       printf "prepare_record_for_matching 4e\n%!"; *)
(*       let product_groups = CommonSubstructure.find_isomorphic_components graph (Collection.of_list (Xlist.map products (fun product -> Collection.of_list (IntSet.to_list product.ids)))) empty_labels in *)
(*       printf "prepare_record_for_matching 4f\n%!"; *)
(**      let groups = Collection.to_list reactant_groups @ Collection.to_list product_groups in
      let _ = Xlist.fold groups 1 (fun n g ->
        Int.fold 1 (Collection.size g) n (fun n i ->
          if n*i > 1000000 then failwith "too many candidates creating during preparing for matching" else
          n*i)) in
      let component_perms = Collection.generate_product_permutations (Collection.of_list (Xlist.rev_map groups (fun c -> c,c))) in**)
(*       printf "prepare_record_for_matching 5\n%!"; *)
    {empty_reaction with
     reactants = reactants;
     products = products;
     ids_reactants = ids_reactants;
     ids_products = ids_products;
     reactant_ids = reactant_ids;
     product_ids = product_ids;
     reactant_groups = Collection.empty(*reactant_groups*);
     product_groups = Collection.empty(*product_groups*);
(*      component_perms = (*component_perms*)Collection.empty; *)
     one_atom_reactants = IntSet.empty(*remove_hydrogens graph one_atom_reactants*);
     one_atom_products = IntSet.empty(*remove_hydrogens graph one_atom_products*);
     cycles = IntSet.empty(*cycles*);
     graph = graph;
     original_graph = graph;
     broken_bonds = 0;
     unbreakable = Collection.empty;
     empty_labels = empty_labels;
     reaction_size = ids_products;
     reactant_stoi = IntSet.empty(*reactant_stoi*);
     product_stoi = IntSet.empty(*product_stoi*);
     msg = [];
     record = r}, messages

let prepare_record_for_matching_no_stoi r messages =
   (* let r = find_aux_solvents r in  *)
(*   printf "prepare_record_for_matching 1\n%!"; *)
  let reactant_smiles,product_smiles = parse_smile_reaction r.reaction_smile in
  let reactants,products = Pair.list_map (reactant_smiles,product_smiles) (make_molecule false) in
  let reactants,ids_reactants = assign_unique_ids_list 1 reactants in
  let products,ids_products = assign_unique_ids_list ids_reactants products in
  let graph = make_atom_graph ids_products (Xlist.map (reactants @ products) (fun m -> m.tree)) in
  (* let graph = collapse_hydrogens_and_fluors graph in
  let reactants = Xlist.map reactants (fun m -> let ids,hf_ids = select_non_hf_atoms graph m.ids in {m with ids=ids; hf_ids=hf_ids}) in
  let products = Xlist.map products (fun m -> let ids,hf_ids = select_non_hf_atoms graph m.ids in {m with ids=ids; hf_ids=hf_ids}) in *)
  let reactant_ids = Xlist.fold reactants IntSet.empty (fun set m -> IntSet.union set m.ids) in
  let product_ids = Xlist.fold products IntSet.empty (fun set m -> IntSet.union set m.ids) in
(*       printf "prepare_record_for_matching 4b\n%!"; *)
(*    let (reactant_ids,one_atom_reactants),(product_ids,one_atom_products) = Pair.list_fold (reactant_trees,product_trees) ((reactant_ids,IntSet.empty),(product_ids,IntSet.empty)) (fun (setm,seto) t ->
      let id = is_one_atom_tree t in
      if id = -1 then setm, seto else IntSet.remove setm id, IntSet.add seto id) in*)
(*  let cycles = Xlist.fold (reactants @ products) IntSet.empty (fun cycles m ->
    Xlist.fold m.cycles cycles (fun cycles (_,set) -> IntSet.union set cycles)) in*)
(*       printf "prepare_record_for_matching 4c\n%!"; *)
  let empty_labels = Labels.make ids_products in
(*       printf "prepare_record_for_matching 4d\n%!"; *)
(*       let reactant_groups = CommonSubstructure.find_isomorphic_components graph (Collection.of_list (Xlist.map reactants (fun reactant -> Collection.of_list (IntSet.to_list reactant.ids)))) empty_labels in *)
(*       printf "prepare_record_for_matching 4e\n%!"; *)
(*       let product_groups = CommonSubstructure.find_isomorphic_components graph (Collection.of_list (Xlist.map products (fun product -> Collection.of_list (IntSet.to_list product.ids)))) empty_labels in *)
(*       printf "prepare_record_for_matching 4f\n%!"; *)
(**      let groups = Collection.to_list reactant_groups @ Collection.to_list product_groups in
      let _ = Xlist.fold groups 1 (fun n g ->
        Int.fold 1 (Collection.size g) n (fun n i ->
          if n*i > 1000000 then failwith "too many candidates creating during preparing for matching" else
          n*i)) in
      let component_perms = Collection.generate_product_permutations (Collection.of_list (Xlist.rev_map groups (fun c -> c,c))) in**)
(*       printf "prepare_record_for_matching 5\n%!"; *)
    {empty_reaction with
     reactants = reactants;
     products = products;
     ids_reactants = ids_reactants;
     ids_products = ids_products;
     reactant_ids = reactant_ids;
     product_ids = product_ids;
     reactant_groups = Collection.empty(*reactant_groups*);
     product_groups = Collection.empty(*product_groups*);
(*      component_perms = (*component_perms*)Collection.empty; *)
     one_atom_reactants = IntSet.empty(*remove_hydrogens graph one_atom_reactants*);
     one_atom_products = IntSet.empty(*remove_hydrogens graph one_atom_products*);
     cycles = IntSet.empty(*cycles*);
     graph = graph;
     original_graph = graph;
     broken_bonds = 0;
     unbreakable = Collection.empty;
     empty_labels = empty_labels;
     reaction_size = ids_products;
     reactant_stoi = IntSet.empty(*reactant_stoi*);
     product_stoi = IntSet.empty(*product_stoi*);
     msg = [];
     record = r}, messages

let count_reactants_atoms r =
  let reactant_smiles,product_smiles = parse_smile_reaction r.reaction_smile in
  Xlist.fold reactant_smiles 0 (fun n reactant ->
    let m = make_molecule false reactant in
	let qmap = calculate_stoi2 StringQMap.empty m.tree in
	StringQMap.fold qmap n (fun n name v -> if name = "H" then n else n+v))
     
let assign_obligatory obligatory molecules =
  let molecules,_ = Xlist.fold molecules ([],1) (fun (molecules,n) m ->
    (if Xlist.mem obligatory n then {m with obligatory=true} else m) :: molecules, n+1) in
  List.rev molecules

let prepare_pattern (name,smiles,obligatory_reactants,obligatory_products) =
  let reactant_smiles,product_smiles = parse_smile_reaction smiles in
  let reactants,products = Pair.list_map (reactant_smiles,product_smiles) make_molecule_pattern in
  let reactants = assign_obligatory obligatory_reactants reactants in
  let products = assign_obligatory obligatory_products products in
  let pattern_size = 1 + (Xlist.fold reactants 0 (fun n m -> max n (IntSet.max_elt m.ids))) in
  let reactant_graph = make_atom_graph pattern_size (Xlist.map reactants (fun m -> m.tree)) in
  let product_graph = make_atom_graph pattern_size (Xlist.map products (fun m -> m.tree)) in
  let empty_labels = Labels.make pattern_size in
  {empty_reaction with
    reactants = reactants;
    products = products;
    ids_reactants = pattern_size;
    reactant_graph = reactant_graph;
    product_graph = product_graph;
    empty_labels = empty_labels;
    pat_name = name;
    smarts = smiles}


(****
let int_of_valence = function
    Alifatic -> 0
  | Aromatic -> 1



let rec latex_chemfig_reid_of_node labels angles reid = function
    b,SAtom(p,l) ->
(*        let id = try IntMap.find reid p.sid with Not_found -> p.sid in *)
       Printf.sprintf "%s[:%d]{\\textcolor{%s}{%s%s}}%s"
                    (latex_chemfig_of_bond b)
                    angles.(p.sid)
                    (latex_color (Labels.get labels p.sid))
                    p.sname (if p.sid = 0 then "" else "_{" ^ string_of_int (try IntMap.find reid p.sid with Not_found -> p.sid) ^ "}")
                    (latex_chemfig_reid_of_list labels angles reid l)
  | b,SLink s -> Printf.sprintf "?[%s,%d]" s (int_of_bond b)

and latex_chemfig_reid_of_list labels angles reid = function
    [] -> ""
  | [x] -> latex_chemfig_reid_of_node labels angles reid x
  | x :: l -> "(" ^ latex_chemfig_reid_of_node labels angles reid x ^ ")" ^ latex_chemfig_reid_of_list labels angles reid l

let latex_chemfig_reid_of labels angles reid tree =
  "\\chemfig{" ^
  (match tree with
    SAtom(p,l) ->
       let id = try IntMap.find reid p.sid with Not_found -> p.sid in
       Printf.sprintf "{\\textcolor{%s}{%s%s}}%s"
                    (latex_color (Labels.get labels p.sid))
                    p.sname (if id = 0 then "" else "_{" ^ string_of_int id ^ "}")
                    (latex_chemfig_reid_of_list labels angles reid l)
  | SLink _ -> failwith "latex_chemfig_of")
  ^ "}"

let rec latex_chemfig_reid2_of_node id labels center_bonds center_links_bonds angles reid = function
    b,SAtom(p,l) ->
(*        let id = try IntMap.find reid p.sid with Not_found -> p.sid in *)
       Printf.sprintf "%s[:%d%s]{\\textcolor{%s}{%s%s}}%s"
                    (latex_chemfig_of_bond b)
                    angles.(p.sid)
                    (try if IntSet.mem (IntMap.find center_bonds (min p.sid id)) (max p.sid id) then ",,,,red" else "" with Not_found -> "")
                    (latex_color (Labels.get labels p.sid))
                    p.sname (if p.sid = 0 then "" else "_{" ^ string_of_int (try IntMap.find reid p.sid with Not_found -> p.sid) ^ "}")
                    (latex_chemfig_reid2_of_list p.sid labels center_bonds center_links_bonds angles reid l)
  | b,SLink s -> Printf.sprintf "?[%s,%d%s]" s (int_of_bond b) (if StringSet.mem center_links_bonds s then ",red" else "")

and latex_chemfig_reid2_of_list id labels center_bonds center_links_bonds angles reid = function
    [] -> ""
  | [x] -> latex_chemfig_reid2_of_node id labels center_bonds center_links_bonds angles reid x
  | x :: l -> "(" ^ latex_chemfig_reid2_of_node id labels center_bonds center_links_bonds angles reid x ^ ")" ^ latex_chemfig_reid2_of_list id labels center_bonds center_links_bonds angles reid l

let latex_chemfig_reid2_of labels center_bonds angles reid tree =
  "\\chemfig{" ^
  (match tree with
    SAtom(p,l) ->
       let id = try IntMap.find reid p.sid with Not_found -> p.sid in
       let links = Xlist.fold l StringMap.empty (fun links a -> make_links_map links p a) in
       let center_links_bonds = StringMap.fold links StringSet.empty (fun center_links_bonds t -> function
            (_,p,[_,q]) -> (try if IntSet.mem (IntMap.find center_bonds (min p.sid q.sid)) (max p.sid q.sid) then StringSet.add center_links_bonds t else center_links_bonds with Not_found -> center_links_bonds)
          | _ -> failwith "latex_chemfig_reid2_of 3") in
       Printf.sprintf "{\\textcolor{%s}{%s%s}}%s"
                    (latex_color (Labels.get labels p.sid))
                    p.sname (if id = 0 then "" else "_{" ^ string_of_int id ^ "}")
                    (latex_chemfig_reid2_of_list p.sid labels center_bonds center_links_bonds angles reid l)
  | SLink _ -> failwith "latex_chemfig_of")
  ^ "}"

let rec change_ids map = function
  | SAtom(p,q) ->
     let p = {p with id=try IntMap.find map p.sid with Not_found -> 0} in
     SAtom(p,Xlist.map q (fun (b,a) -> b,change_ids map a))
  | SLink t -> SLink t

(*let rec change_tree_labels ids = function
  | SAtom(p,q) ->
     let p = {p with label=try IntMap.find ids p.sid with Not_found -> 0} in
     let q = Xlist.map q (fun (b,a) -> b, change_tree_labels ids a) in
     SAtom(p,q)
  | SLink t -> SLink t*)

let rec count_tree_size = function
  | SAtom(_,q) -> Xlist.fold q 1 (fun n (_,a) -> n + count_tree_size a)
  | SLink _ -> 0

(*let rec count_unlabelled_atoms n = function
  | SAtom(p,q) ->
(*       Printf.printf "%d %s\n" n (string_of_tree_labels (SAtom(p,q))); *)
      Xlist.fold q (if p.label=0 then n+1 else n) (fun n (_,a) -> count_unlabelled_atoms n a)
  | SLink _ -> n

let rec get_unlabelled_atoms set = function
  | SAtom(p,q) ->
(*       Printf.printf "%d %s\n" n (string_of_tree_labels (SAtom(p,q))); *)
      Xlist.fold q (if p.label=0 then IntSet.add set p.sid else set) (fun set (_,a) -> get_unlabelled_atoms set a)
  | SLink _ -> set*)

let rec rename_atoms names = function
  | SAtom(p,q) ->
     let p = if IntMap.mem names p.sid then
        {name = IntMap.find names p.sid; valence=Alifatic; charge=0; chirality=Unspecified; hydrogens=(-1); id=p.sid}
       else p in
     let q = Xlist.map q (fun (b,a) -> b, rename_atoms names a) in
     SAtom(p,q)
  | SLink t -> SLink t

let rec rename_atoms2 name new_name = function
  | SAtom(p,q) ->
     let p = if p.sname = name then
        {p with name = new_name}
       else p in
     let q = Xlist.map q (fun (b,a) -> b, rename_atoms2 name new_name a) in
     SAtom(p,q)
  | SLink t -> SLink t

let rec merge_atoms names = function
  | SAtom(p,q) ->
     let q = Xlist.map q (fun (b,a) -> b, merge_atoms names a) in
     let q = List.rev (Xlist.fold q [] (fun q -> function
         b,SAtom(p2,q2) -> if p.sname = p2.sname && StringSet.mem names p.sname then q2 @ q else (b,SAtom(p2,q2)) :: q
       | b,SLink t -> (b,SLink t) :: q)) in
     SAtom(p,q)
  | SLink t -> SLink t

let rec remove_links names = function
  | SAtom(p,q) ->
     let q = Xlist.map q (fun (b,a) -> b, remove_links names a) in
     let links = Xlist.fold q StringQMap.empty (fun links -> function
         _,SAtom _ -> links
       | _,SLink t -> StringQMap.add links t) in
     let q = List.rev (Xlist.fold q [] (fun q -> function
         b,SAtom(p2,q2) -> (b,SAtom(p2,q2)) :: q
       | b,SLink t -> if StringQMap.find links t = 2 then q else (b,SLink t) :: q)) in
     SAtom(p,q)
  | SLink t -> SLink t

let rec change_root_node name root = function
    SAtom(p,l) ->
      if p.sname = name then [SAtom(p,root @ l)] else
      change_root_list name root p [] l
  | SLink _ -> []

and change_root_list name root p rev = function
    (b,a) :: l -> (change_root_node name [b,SAtom(p,root @ List.rev rev @ l)] a) @ (change_root_list name root p ((b,a) :: rev) l)
  | [] -> []

let change_root name t =
  match change_root_node name [] t with
    [t] -> t
  | _ -> (*failwith "change_root"*)t

let rec label_change labels label rev = function
    [] -> raise Not_found
  | (b,SLink s) :: l -> label_change labels label ((b,SLink s) :: rev) l
  | (b,SAtom(p,l2)) :: l ->
       let label2 = try Labels.find labels p.sid with Not_found -> failwith "label_change" in
       if label = label2 then label_change labels label ((b,SAtom(p,l2)) :: rev) l
       else List.rev rev @ l, SAtom(p,l2), label2, p.sid

let rec split_labeled_node labels root = function
    SAtom(p,l) ->
      (try
        let label = try Labels.find labels p.sid with Not_found -> failwith "split_labeled_node" in
        let l,a,label2,id2 = label_change labels label [] l in
        [SAtom({name="X"; valence=Alifatic; charge=0; chirality=Unspecified; hydrogens=(-1); id=0},[Default,SAtom(p,root @ l)]), label, p.sid,
        SAtom({name="X"; valence=Alifatic; charge=0; chirality=Unspecified; hydrogens=(-1); id=0},[Default,a]), label2, id2]
      with Not_found -> split_labeled_list labels root p [] l)
  | SLink s -> []

and split_labeled_list labels root p rev = function
    (b,a) :: l -> (split_labeled_node labels [b,SAtom(p,root @ List.rev rev @ l)] a) @ (split_labeled_list labels root p ((b,a) :: rev) l)
  | [] -> []

let split_labeled labels t =
  match split_labeled_node labels [] t with
    [t] -> t
  | _ -> failwith "split_labeled"

let check_smile_correctness_list l =
  Xlist.iter l (fun s ->
    match parse_smile_molecule "" s with
      SLink s -> failwith s
    | _ -> ())

(*let check_smile_correctness () =
  Import.fold_rxn () (fun () r ->
    try
      check_smile_correctness_list r.reactant_smiles;
      check_smile_correctness_list r.product_smiles;
      check_smile_correctness_list r.solvent_smiles
    with Failure s -> (
      Printf.printf "%s\n" s;
      print_endline (Import.string_of_record_fmt r)))  *)


(* sprawdzenie poprawności smilesów *)
(* let _ = check_smile_correctness () *)
(* wynik w pliku smile_correctness.txt *)

(*let check_smiles_translation () =
  let molecules = Import.load_molecules () in
  Xlist.iter molecules (fun smile ->
    let tree = parse_smile_molecule "" smile in
    let smile2 = string_of_tree_std_id tree in
    if smile <> smile2 then Printf.printf "%s --> %s\n%!" smile smile2)*)


(* sprawdzenie poprawności konwersji pomiędzy smiles i smiles_tree *)
(* let _ = check_smiles_translation () *)
****)
