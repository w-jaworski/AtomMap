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
open Printf


(*let calculate_stoi_reaction r =
  let reactant_smiles,product_smiles = Smiles.parse_smile_reaction r.reaction_smile in
  let reactant_trees,product_trees = Pair.list_map (reactant_smiles,product_smiles) (fun s -> (*List.hd*) (Smiles.manage_aromaticity (Smiles.parse_smile_molecule r.rxn_id s))) in
  let _,reactant_trees,size = Smiles.assign_unique_ids_and_translate_to_tree_list 0 reactant_trees in
  let _,product_trees,size = Smiles.assign_unique_ids_and_translate_to_tree_list size product_trees in
  calculate_stoi reactant_trees product_trees*)



(*let parse_reaction r =
  let reactant_smiles,product_smiles = Smiles.parse_smile_reaction r.reaction_smile in
  {r with reactant_smiles = reactant_smiles;
          product_smiles = product_smiles}*)

(*let rec remove_aux_solvents rev s = function
    [] -> failwith "remove_aux_solvents"
  | x :: l -> if x = s then List.rev rev @ l else remove_aux_solvents (x :: rev) s l

let find_aux_solvents r =
  let l = Pair.string_bucket_group (Pair.singleton r.reactant_smiles r.product_smiles) (fun x -> x) (fun x -> x) in
  let aux_solvents = Pair.fold l [] (fun l (x,y) ->
    if Xlist.size x < Xlist.size y then x @ l else y @ l) in
  let reactants,products = Pair.list_fold (aux_solvents,aux_solvents) (r.reactant_smiles,r.product_smiles) (fun l s -> remove_aux_solvents [] s l) in
  {r with aux_solvent_smiles = aux_solvents;
          reactant_smiles = reactants;
          product_smiles = products}*)
(****

let remove_hydrogens graph set =
  IntSet.fold set IntSet.empty (fun set i ->
    if (fst graph.(i)).name = "H" then set else IntSet.add set i)

let remove_hydrogens2 graph =
  Array.map (fun (p,l) ->
    if p.name = "H" then p,[] else
    let l = List.rev (Xlist.fold l [] (fun l (b,q) -> if q.name = "H" then l else (b,q) :: l)) in
    p,l) graph

let rec is_one_atom_tree = function
    SAtom(p,[]) ->(* Printf.printf "p.id=%d\n%!" p.id;*) p.id
  | SAtom(p,l) -> (*print_endline (Smiles.string_of_smile_tree (SAtom(p,l)));*)
             if p.name = "H" && Xlist.size l = 1 then is_one_atom_tree (snd (List.hd l)) else
             let b = Xlist.fold l true (fun b -> function
                _,SAtom(q,[]) -> if q.name = "H" then b else false
              | _ -> false) in
             (*let n =*) if b then p.id else -1 (*in
             Printf.printf "nn=%d p.id=%d\n%!" n p.id; n*)
  | _ -> -1

let rec remove_one_atom_trees = function
    [] -> []
  | t :: l -> if is_one_atom_tree t = -1 then t :: (remove_one_atom_trees l) else remove_one_atom_trees l

let no_reactants r =
  Xlist.size (remove_one_atom_trees r.reactant_trees)

let no_products r =
  Xlist.size (remove_one_atom_trees r.product_trees)
****)

(*let disambiguate_molecule empty_labels subtasks =
  let subtasks_ann = Collection.annotate subtasks (fun r ->
    r.graph, Hash.make_hash Hash.full_hash_key 4 r.graph r.empty_labels) in
  let label_buckets =
    try Collection.quotient subtasks_ann (fun x y -> Isomorphism.are_isomorphic_alt_lists true 4 x y component_perms)
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during partitioning") in
(*  ignore (Xlist.fold (Collection.to_list label_buckets) 1 (fun i l ->
    Xlist.iter (Collection.to_list l) (fun (r,(_,h)) ->
      printf "d%d unbreakable=%s\n%!" i (String.concat " " (Xlist.map r.unbreakable (fun (i,j) -> sprintf "%d-%d" i j)));
      print_graph r.graph h);
    i+1));*)
  Collection.to_list (Collection.deannotate (Collection.reverse_projection label_buckets))*)



(**********************************************************************************************)

(****
let rec iter_label_reaction min_level r labels = function
    0 -> [labels]
  | 1 ->
        let l = Division.label_reaction min_level r labels in
        Xlist.iter l (fun labels ->
          Printf.printf "%d reactants %s\n" min_level (String.concat " " (Xlist.map r.reactant_trees (Smiles.string_of_tree_terminal_colors labels)));
          Printf.printf "%d products  %s\n%!\n" min_level (String.concat " " (Xlist.map r.product_trees (Smiles.string_of_tree_terminal_colors labels))));
        l
  | n ->
     (try
        let l = Division.label_reaction (min_level+n-1) r labels in
        Xlist.iter l (fun labels ->
          Printf.printf "%d reactants %s\n" (min_level+n-1) (String.concat " " (Xlist.map r.reactant_trees (Smiles.string_of_tree_terminal_colors labels)));
          Printf.printf "%d products  %s\n%!\n" (min_level+n-1) (String.concat " " (Xlist.map r.product_trees (Smiles.string_of_tree_terminal_colors labels))));
        List.flatten (Xlist.map l (fun labels -> iter_label_reaction min_level r labels (n-1)))
      with Failure "no exclusion schemes found" -> iter_label_reaction min_level r labels (n-1))

let get_level = function
    [] -> Parse 0
  | P 0 :: _ -> Parsed 0
  | P n :: _ -> Parse (n-1)
  | N _ :: P n :: _ -> Parsed n
  | A _ :: P n :: _ -> ParsedAmb n
  | N 4 :: _ -> FailNone
  | A 4 :: _ -> FailAmb
  | N n :: _ -> Parse (n+1)
  | A n :: _ -> Parse (n+1)
  | _ -> failwith "get_level"
(*   | _ -> Unknown *)

let get_level2 = function
    [] -> Parse 4
  | P 1 :: _ -> Parsed 1
  | P n :: _ -> BorderParse (n-1)
  | B n :: _ -> Parse n
  | N 4 :: _ -> Parse 3
  | N 1 :: _ -> FailNone
  | N n :: _ -> BorderParse (n-1)
  | A _ :: _ -> FailAmb
  | NB n :: _ -> Parse n
  | AB _ :: _ -> FailAmb
(*   | _ -> Unknown *)
****)
(****
let rec multilevel_label_reaction r history labels_list =
  match get_level history with
    Parse level ->
      (try
        let alt_labels_list =
          if labels_list = [] then Division.label_reaction level r r.empty_labels
          else List.flatten (Xlist.map labels_list (fun labels ->
            try Division.label_reaction level r labels with Division.SolutionNotFound _ -> [])) in
        if alt_labels_list = [] then raise (Division.SolutionNotFound "cummulative") else
        multilevel_label_reaction r ((P level) :: history) alt_labels_list
      with Division.Ambiguity _ -> multilevel_label_reaction r ((A level) :: history) labels_list
         | Division.SolutionNotFound _ -> multilevel_label_reaction r ((N level) :: history) labels_list)
  | result -> result, List.rev history, labels_list

let rec multilevel_label_reaction2 r history labels_list =
  match get_level2 history with
    Parse level ->
      (try
        let alt_labels_list =
          if labels_list = [] then Division.label_reaction level r r.empty_labels
          else List.flatten (Xlist.map labels_list (fun labels ->
            try Division.label_reaction level r labels with Division.SolutionNotFound _ -> [])) in
        if alt_labels_list = [] then raise (Division.SolutionNotFound "cummulative") else
        multilevel_label_reaction2 r ((P level) :: history) alt_labels_list
      with Division.Ambiguity _ -> multilevel_label_reaction2 r ((A level) :: history) labels_list
         | Division.SolutionNotFound _ -> multilevel_label_reaction2 r ((N level) :: history) labels_list)
  | BorderParse level ->
      (try
        let alt_labels_list =
          if labels_list = [] then Division.label_border_reaction level r r.empty_labels
          else List.flatten (Xlist.map labels_list (fun labels ->
            try Division.label_border_reaction level r labels with Division.SolutionNotFound _ -> [])) in
        if alt_labels_list = [] then raise (Division.SolutionNotFound "cummulative") else
        multilevel_label_reaction2 r ((B level) :: history) alt_labels_list
      with Division.Ambiguity _ -> multilevel_label_reaction2 r ((AB level) :: history) labels_list
         | Division.SolutionNotFound _ -> multilevel_label_reaction2 r ((NB level) :: history) labels_list)
  | result -> result, List.rev history, labels_list
****)
(****
let string_of_result = function
    Parse level -> "Parse_" ^ string_of_int level
  | BorderParse level -> "BorderParse_" ^ string_of_int level
  | Parsed level -> "Parsed_" ^ string_of_int level
  | ParsedAmb level -> "ParsedAmb_" ^ string_of_int level
  | FailAmb -> "FailAmb"
  | FailNone -> "FailNone"
  | Unknown -> "Unknown"

(**********************************************************************************************)

let get_border_schemas ll (labels,division) =
  let divisions = Pair.get_border_schemas division ll in
  Xlist.rev_map divisions (fun division -> Labels.add_list_of_list_pairs labels division, division)

let rec match_border_atoms r (labels,matching) =
  let hash = Hash.make_border_hash r.graph labels in
  let hash_buckets = Hash.hash_bucket_group hash r.reactant_ids r.product_ids in
  try
    let no_ex_sch = Pair.int_log (Pair.number_of_exclusion_schemes hash_buckets) in
    if no_ex_sch > 3 then failwith "too many exclusion schemes" else
    let alt_labels_list = get_border_schemas hash_buckets (labels,matching) in
    Xlist.fold alt_labels_list [] (fun alt_labels_list (labels,matching) ->
      alt_labels_list @ match_border_atoms r (labels,matching))
  with Pair.Empty -> [labels,matching]

let match_border_atoms2 r (labels,matching) =
  let hash = Hash.make_border_hash r.graph labels in
  let hash_buckets = Hash.hash_bucket_group hash r.reactant_ids r.product_ids in
  try
    let no_ex_sch = Pair.int_log (Pair.number_of_exclusion_schemes hash_buckets) in
    if no_ex_sch > 3 then failwith "too many exclusion schemes" else
    let alt_labels_list = get_border_schemas hash_buckets (labels,matching) in
    alt_labels_list
  with Pair.Empty -> [labels,matching]

let match_border_atoms3 r labels =
  let hash = Hash.make_border_hash r.graph labels in
  let hash_buckets = Hash.hash_bucket_group hash r.reactant_ids r.product_ids in
  let b = Pair.fold hash_buckets true (fun b (l1,l2) -> b && (Xlist.size l1 = Xlist.size l2)) in
  if b then
    get_border_schemas hash_buckets (labels,Pair.empty)
  else []

let disambiguate_matchings hash matchings =
  let map = Xlist.fold matchings IntMap.empty (fun map matching ->
    let l = Xlist.fold matching [] (fun l (v,w) ->
      if hash.(v) <> hash.(w) then v :: w :: l else l) in
    IntMap.add_inc map (Xlist.size l) [matching(*,l*)] (fun ll -> (matching(*,l*)) :: ll)) in
  IntMap.fold map ([],1000000) (fun (x,v) k y -> if v <= k then x,v else y,k)

let match_reaction r labels division =
  let hash = Hash.make_hash Hash.hash_key 1 r.graph labels in
  let unlabelled_reactants = AtomGraph.get_unlabelled_atoms r.reactant_ids labels in
  let unlabelled_products = AtomGraph.get_unlabelled_atoms r.product_ids labels in
  let selected_reactants = IntSet.fold r.reactant_ids IntSet.empty (fun set i -> if Labels.mem labels i then IntSet.add set i else set) in
  let selected_products = IntSet.fold r.product_ids IntSet.empty (fun set i -> if Labels.mem labels i then IntSet.add set i else set) in
  let matchings : (int * int) list list = Division.extend_matching3 hash.(1) r.graph selected_reactants selected_products unlabelled_reactants unlabelled_products division in
  if matchings = [] then (*(print_endline ("empty matchings " ^ r.record.rxn_id);*) [] else
  let hash = Hash.make_hash Hash.full_hash_key 1 r.graph r.empty_labels in
  let (matchings2 : (int * int) list list),k = disambiguate_matchings hash.(1) matchings in
(*           Printf.printf "A id=%s %d %d k=%d\n" r.record.rxn_id (Xlist.size matchings) (Xlist.size matchings2) k; *)
  [matchings2,k]

****)
