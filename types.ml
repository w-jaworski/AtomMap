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

open Xstd

(* TODO:
  konwersja experimental_yield, min_temp, max_temp
  pomijam chiralność i slashe w modelu grafowym
*)

type bond = Default | Single | Slash | Backslash | Double | Triple | Aromatic

type valence = Alifatic | Aromatic

type chirality = Unspecified | Clockwise | AntiClockwise

type smile_atom_props = {sname: string;
                   svalence: valence;
                   scharge: int;
                   schirality: chirality;
                   shydrogens: int;
                   sid: int}

type atom_props = {name: string;
(*                    valence: valence;  *)
                   charge: int;
                   chirality: chirality;
(*                    hydrogens: int;  *)
                   id: int;
                   hydrogens: int list;
                   fluors: int list;
                   (*label: int*)}

type smile_tree =
    SAtom of smile_atom_props * (bond * smile_tree) list
  | SLink of string

type atom_tree =
    Atom of atom_props * (bond * atom_tree) list
  | Link of string

(* type molecule_type = Aux of int | Half of int *)

type graph = (atom_props * (bond * atom_props) list) array

type hash = string array

type record = {rxn_id: string; patent_id: string; path: string;
               reactant_smiles: string list;
               product_smiles: string list;
               solvent_smiles: string list;
               solvents: string;
               experimental_yield: string; experimental_yield_range: string; experimental_yield_calc: string;
               temp: string; min_temp: string; max_temp: string;
(*               missing_reactant_smiles: (molecule_type * string) list;
               missing_product_smiles: (molecule_type * string) list;*)
               reaction_smile: string;
               added_reactant_smiles: string list;
               added_product_smiles: string list;
               bibref: string;
               is_correct: bool;
               (*aux_reactant_smiles: string list;
               aux_product_smiles: string list;
               aux_solvent_smiles: string list*)}

type molecule = {smiles: string;
                 smile_tree: smile_tree;
                 tree: atom_tree;
                 ids: IntSet.t;
                 hf_ids: IntSet.t;
                 stoi: bool;
                 obligatory: bool;
                 cycles: (Xstd.StringSet.t * Xstd.IntSet.t) list}

type reaction = {reactants: molecule list;
                 products: molecule list;
                 ids_reactants: int;
                 ids_products: int;
                 reactant_ids: IntSet.t;
                 product_ids: IntSet.t;
                 reactant_stoi: IntSet.t;
                 product_stoi: IntSet.t;
                 reactant_groups: int Collection.collection Collection.collection Collection.collection;
                 product_groups: int Collection.collection Collection.collection Collection.collection;
(*                  component_perms: (int Collection.collection * int Collection.collection) Collection.collection Collection.collection; *)
                 one_atom_reactants: IntSet.t;
                 one_atom_products: IntSet.t;
                 cycles: IntSet.t;
                 original_graph: graph;
                 graph: graph;
                 broken_bonds: int;
                 unbreakable: (int * int) Collection.collection;
                 empty_labels: Labels.t;
                 reaction_size: int;
                 msg: string list;
                 record: record;
                 reactant_graph: graph;
                 product_graph: graph;
                 pat_name: string;
                 smarts: string}

let empty_record =   {rxn_id = ""; patent_id = ""; path = "";
                      reactant_smiles = [];
                      product_smiles = [];
                      solvent_smiles = [];
                      solvents = "";
(*                      missing_reactant_smiles = [];
                      missing_product_smiles = [];*)
(*                      aux_reactant_smiles = [];
                      aux_product_smiles = [];
                      aux_solvent_smiles = []; *)
                      reaction_smile = "";
                      added_reactant_smiles = [];
                      added_product_smiles = [];
                      experimental_yield = "";  experimental_yield_range = ""; experimental_yield_calc = "";
                      temp = ""; min_temp = ""; max_temp = "";
                      is_correct = false;
                      bibref = ""}

let empty_reaction = {
     reactants = [];
     products = [];
     ids_reactants = 0;
     ids_products = 0;
     reactant_ids = IntSet.empty;
     product_ids = IntSet.empty;
     reactant_groups = Collection.empty;
     product_groups = Collection.empty;
     one_atom_reactants = IntSet.empty;
     one_atom_products = IntSet.empty;
     cycles = IntSet.empty;
     graph = [| |];
     original_graph = [| |];
     broken_bonds = 0;
     unbreakable = Collection.empty;
     empty_labels = Labels.make 1;
     reaction_size = 0;
     reactant_stoi = IntSet.empty;
     product_stoi = IntSet.empty;
     msg = [];
     record = empty_record;
     reactant_graph = [| |];
     product_graph = [| |];
     pat_name = "";
     smarts = ""}


type pattern_atom_props = {names: string list;
                   pid: int;
                  }

type pattern_graph = (pattern_atom_props * (bond * pattern_atom_props) list) array

type smart_tree =
    PAtom of pattern_atom_props * (bond * smart_tree) list
  | PLink of int

type pattern = {smart: string;
                smart_tree: smart_tree;
                }

type rearrangement_pattern = {
  rname: string;
  patterns: pattern list;
  pattern_graph: pattern_graph}

(****type node =
    Atom of int *  string * smile_atom_props (* id * name *)
  | Complex of node list * edge list

and path =
    Path of edge * node * path list
  | Link of edge * string
  | Loop of edge
  | Leaf of edge * node * path list
  | Root of node * path list
  | Ion of string

and edge =
    Bond of int * bond * int
  | Edge of (edge * node * path list) list * edge
  | Parallel of edge list****)

type message_from_overseer =
    Work_with of int * (bool * record)
  | Kill_yourself

type message_to_overseer =
    Ready_to_work of string
  | Work_done of
      string  * (int * record * string list * (reaction * Labels.t * int IntMap.t IntMap.t * int IntMap.t IntMap.t * IntMap.key IntMap.t) list)



let print_intqmap qmap =
  let sum = IntQMap.fold qmap 0 (fun n _ v -> n+v) in
  Printf.printf " key quantity percentage sum_percent\n";
  let _ = IntQMap.fold qmap 0 (fun n k v ->
    Printf.printf "%4d %8d %.8f %.8f\n" k v (float v /. float sum) (float (n+v) /. float sum);
    n+v) in
  Printf.printf "total quantity: %d\n" sum

let print_stringqmap qmap =
  let sum = StringQMap.fold qmap 0 (fun n _ v -> n+v) in
  StringQMap.iter qmap (fun k v ->
    Printf.printf "%8d %.8f %s\n" v (float v /. float sum) k);
  Printf.printf "total quantity: %d\n" sum

let string_of_int_list l =
  Printf.sprintf "[%s]" (String.concat ";" (Xlist.map l string_of_int))

let string_of_selection l =
  "[" ^ String.concat ";" (Xlist.map (List.sort compare l) string_of_int) ^ "]"

let string_of_selection_list l =
  "[" ^ String.concat ";" (Xlist.map l string_of_selection) ^ "]"


exception Timeout
let time = ref (Sys.time ())
let timeout = 1200.

let logfile = open_out_gen [Open_wronly; Open_append; Open_creat] ((6*8+4)*8+4) "results/worker.log"
