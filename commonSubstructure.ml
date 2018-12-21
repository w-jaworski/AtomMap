(*
 *  AtomMap: maps atoms in chemical reactions
 *  Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl> 
 *  Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw                

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
(****
type t = int Pair.t
****)
let rec get_connected_component removed graph found v =
  if IntSet.mem removed v then found,removed else
  Xlist.fold (snd graph.(v)) (v :: found,IntSet.add removed v) (fun (found,removed) (_,p) ->
    get_connected_component removed graph found p.id)

let get_connected_components e graph selection =
  let removed = IntSet.of_list (Collection.to_list e) in
  let conn = fst (IntSet.fold selection ([],removed) (fun (found,removed) v ->
    let l,removed = get_connected_component removed graph [] v in
    if l = [] then found,removed else
    l :: found, Xlist.fold l removed IntSet.add)) in
  Collection.of_list (Xlist.map conn Collection.of_list)

(****
let extend_matching3 hash g selection1 selection2 removed1 removed2 division =
  let conn1 = get_connected_components removed1 g selection1 in
  let conn2 = get_connected_components removed2 g selection2 in
  Isomorphism.extend_matching3 hash g selection1 selection2 conn1 conn2 division
****)

(* conn - kolekcja spójnych składowych - kolekcja kolekcji wierzchołków
   buckets - kolekcja par kolekcji spójnych składowych
   part - kolekcja kolekcji izomorficznych spójnych składowych
   part_buckets - kolekcja par kolekcji równolicznych kolekcji izomorficznych spójnych składowych
*)

let find_isomorphic_components graph conn labels =
  let hash = Hash.make_hash Hash.hash_key 1 graph labels in
  Collection.quotient conn (fun x y -> Isomorphism.are_isomorphic true 1 hash hash graph graph x y)

let match_isomorphic_connected_components2 e graph selection1 selection2 hash level =
  let conn1 = get_connected_components e graph selection1 in
  let conn2 = get_connected_components e graph selection2 in
  if Collection.size conn1 <> Collection.size conn2 then raise Not_found else
  let buckets = Collection.int_pair_bucket_group conn1 conn2 (fun c -> Collection.size c) in
  if Collection.exists buckets (fun (conn1,conn2) -> Collection.size conn1 <> Collection.size conn2) then raise Not_found else
  let buckets = Collection.flatten_map buckets (fun (conn1,conn2) ->
    Collection.string_pair_bucket_group conn1 conn2 (fun c -> Hash.hash_list_key hash.(level) c)) in
  if Collection.exists buckets (fun (conn1,conn2) -> Collection.size conn1 <> Collection.size conn2) then raise Not_found else
  let part_buckets = Collection.flatten_map buckets (fun (conn1,conn2) ->
    let part1 = Collection.quotient conn1 (fun x y -> Isomorphism.are_isomorphic false level hash hash graph graph x y) in
    let part2 = Collection.quotient conn2 (fun x y -> Isomorphism.are_isomorphic false level hash hash graph graph x y) in
    Collection.int_pair_bucket_group part1 part2 (fun c -> Collection.size c)) in
  if Collection.exists part_buckets (fun (part1,part2) -> Collection.size part1 <> Collection.size part2) then raise Not_found else
  Collection.flatten_map part_buckets (fun (part1,part2) ->
    Collection.combine part1 part2 (fun conn1 conn2 ->
      Isomorphism.are_isomorphic false level hash hash graph graph (Collection.elem conn1) (Collection.elem conn2)))

let match_isomorphic_connected_components excl graph selection1 selection2 hash level =
  Collection.flatten_map excl (fun e ->
    try Collection.singleton (match_isomorphic_connected_components2 e graph selection1 selection2 hash level)
    with Not_found -> Collection.empty)

let match_isomorphic_connected_components_simple graph selection1 selection2 hash level =
  try match_isomorphic_connected_components2 Collection.empty graph selection1 selection2 hash level
  with Not_found -> Collection.empty

exception Ambiguity of string
exception SolutionNotFound of string

let add_label_isomorfic_connected_components labels comp =
  let ll = Pair.of_pair_list (Xlist.map (Pair.to_pair_list comp) (fun (l1,l2) -> List.flatten l1, List.flatten l2)) in
  Labels.add_list_of_list_pairs labels ll
(*  Xlist.fold comp labels_next (fun (labels,next) (l1,l2) ->
    let labels = Xlist.fold l1 labels (fun labels c -> Xlist.fold c labels (fun labels v ->
      if Labels.mem labels v then labels else Labels.set labels v next)) in
    let labels = Xlist.fold l2 labels (fun labels c -> Xlist.fold c labels (fun labels v ->
      if Labels.mem labels v then labels else Labels.set labels v next)) in
    labels, next+1)*)

let label_reaction level r labels =
  let hash = Hash.make_hash Hash.hash_key level r.graph labels in
  let hash_buckets = Hash.hash_bucket_group hash.(level) (Collection.of_list (IntSet.to_list r.reactant_ids)) (Collection.of_list (IntSet.to_list r.product_ids)) in
  let no_ex_sch =
    try Collection.int_log (Collection.number_of_exclusion_schemes hash_buckets)
    with Collection.Empty -> raise (SolutionNotFound "no exclusion schemes found") in
  if no_ex_sch > 3 then raise (Ambiguity "too many exclusion schemes") else
  let ex = Collection.get_exclusion_schemas hash_buckets in
  let ex =
    try match_isomorphic_connected_components ex r.graph r.reactant_ids r.product_ids hash level
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during component matching") in
  if Collection.is_empty ex then raise (SolutionNotFound "no solution found") else
  Xlist.map (Collection.to_list ex) (fun e ->
    let e = Pair.of_pair_list (Xlist.map (Collection.to_list e) (fun (c1,c2) ->
      Xlist.map (Collection.to_list c1) Collection.to_list,
      Xlist.map (Collection.to_list c2) Collection.to_list)) in
    let ll = Pair.of_pair_list (Xlist.map (Pair.to_pair_list e) (fun (l1,l2) -> List.flatten l1, List.flatten l2)) in
    Labels.add_list_of_list_pairs labels ll)

let are_isomorphic_graphs graph1 graph2 selection1 selection2 hash1 hash2 level =
  let conn1 = get_connected_components Collection.empty graph1 selection1 in
  let conn2 = get_connected_components Collection.empty graph2 selection2 in
  if Collection.size conn1 <> Collection.size conn2 then false else
  let buckets = Collection.int_pair_bucket_group conn1 conn2 (fun c -> Collection.size c) in
  if Collection.exists buckets (fun (conn1,conn2) -> Collection.size conn1 <> Collection.size conn2) then false else
  let buckets = Collection.flatten_map buckets (fun (conn1,conn2) ->
    Collection.string_pair_bucket_group2 conn1 conn2 (fun c -> Hash.hash_list_key hash1.(level) c) (fun c -> Hash.hash_list_key hash2.(level) c)) in
  if Collection.exists buckets (fun (conn1,conn2) -> Collection.size conn1 <> Collection.size conn2) then false else
  let part_buckets = Collection.flatten_map buckets (fun (conn1,conn2) ->
    let part1 = Collection.quotient conn1 (fun x y -> Isomorphism.are_isomorphic false level hash1 hash1 graph1 graph1 x y) in
    let part2 = Collection.quotient conn2 (fun x y -> Isomorphism.are_isomorphic false level hash2 hash2 graph2 graph2 x y) in
    Collection.int_pair_bucket_group part1 part2 (fun c -> Collection.size c)) in
  if Collection.exists part_buckets (fun (part1,part2) -> Collection.size part1 <> Collection.size part2) then false else
  let ex = try Collection.flatten_map part_buckets (fun (part1,part2) ->
    Collection.combine part1 part2 (fun conn1 conn2 ->
      Isomorphism.are_isomorphic false level hash1 hash2 graph1 graph2 (Collection.elem conn1) (Collection.elem conn2)))
    with Not_found -> Collection.empty in
  not (Collection.is_empty ex)

let disambiguate_labels level r alt_labels_list =
(*   Printf.printf "disambiguate_labels |alt_labels_list|=%d\n%!" (Xlist.size alt_labels_list); *)
  let alt_labels = Collection.of_list alt_labels_list in
  let alt_labels_ann = Collection.annotate alt_labels (fun labels ->
    r.graph, Hash.make_hash Hash.full_hash_key level r.graph labels) in
  let label_buckets =
    try Collection.quotient alt_labels_ann (fun (x1,(graph1,hash1)) (x2,(graph2,hash2)) ->
      are_isomorphic_graphs graph1 graph2 r.reactant_ids r.reactant_ids hash1 hash2 level &&
      are_isomorphic_graphs graph1 graph2 r.product_ids r.product_ids hash1 hash2 level)
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during partitioning") in
  Collection.to_list (Collection.deannotate (Collection.reverse_projection label_buckets))

(*let disambiguate_labels level r alt_labels_list =
(*   Printf.printf "disambiguate_labels |alt_labels_list|=%d\n%!" (Xlist.size alt_labels_list); *)
  let alt_labels = Collection.of_list alt_labels_list in
  let alt_labels_ann = Collection.annotate alt_labels (fun labels ->
    r.graph, Hash.make_hash Hash.full_hash_key level r.graph labels) in
  let label_buckets =
    try Collection.quotient alt_labels_ann (fun x y -> Isomorphism.are_isomorphic_alt_lists true level x y (Collection.to_list r.component_perms))
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during partitioning") in
  Collection.to_list (Collection.deannotate (Collection.reverse_projection label_buckets))*)

let label_border_reaction level r labels =
  let hash = Hash.make_border_hash2 level r.graph labels in
(*   Int.iter 0 (Array.length hash.(level) - 1) (fun i -> Printf.printf "%d %s\n" i hash.(level).(i)); *)
  let hash_buckets = Hash.hash_bucket_group hash.(level) (Collection.of_list (IntSet.to_list r.reactant_ids)) (Collection.of_list (IntSet.to_list r.product_ids)) in
  let reactant_glob_ex = IntSet.fold r.reactant_ids [] (fun l i -> if hash.(level).(i) = "" && not (Labels.mem labels i) then i :: l else l) in
  let product_glob_ex = IntSet.fold r.product_ids [] (fun l i -> if hash.(level).(i) = "" && not (Labels.mem labels i) then i :: l else l) in
  let glob_ex = Collection.sum (Collection.of_list reactant_glob_ex) (Collection.of_list product_glob_ex) in
(*  StringMap.iter hash_buckets (fun k (l1,l2) ->
    Printf.printf "%s %s %s\n" k (AtomGraph.string_of_int_list l1) (AtomGraph.string_of_int_list l2));*)
  let no_ex_sch =
    try Collection.int_log (Collection.number_of_exclusion_schemes hash_buckets)
    with Collection.Empty -> raise (SolutionNotFound "no exclusion schemes found") in
(*   Printf.printf "no_ex_sch = %d\n" no_ex_sch; *)
  if no_ex_sch > 3 then raise (Ambiguity "too many exclusion schemes") else
  let ex = Collection.get_exclusion_schemas hash_buckets in
(*  Printf.printf "reactant_ex [%s]\n" (String.concat ";" (Xlist.map reactant_ex (fun l -> AtomGraph.string_of_int_list l)));
  Printf.printf "product_ex [%s]\n" (String.concat ";" (Xlist.map product_ex (fun l -> AtomGraph.string_of_int_list l)));*)
  let ex = Collection.map ex (Collection.sum glob_ex) in
(*  Printf.printf "reactant_ex [%s]\n" (String.concat ";" (Xlist.map reactant_ex (fun l -> AtomGraph.string_of_int_list l)));
  Printf.printf "product_ex [%s]\n" (String.concat ";" (Xlist.map product_ex (fun l -> AtomGraph.string_of_int_list l)));*)
  let ex =
    try match_isomorphic_connected_components ex r.graph r.reactant_ids r.product_ids hash level
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during component matching") in
  if Collection.is_empty ex then raise (SolutionNotFound "no solution found") else
  Xlist.map (Collection.to_list ex) (fun e ->
    let e = Pair.of_pair_list (Xlist.map (Collection.to_list e) (fun (c1,c2) ->
      Xlist.map (Collection.to_list c1) Collection.to_list,
      Xlist.map (Collection.to_list c2) Collection.to_list)) in
     let ll = Pair.of_pair_list (Xlist.map (Pair.to_pair_list e) (fun (l1,l2) -> List.flatten l1, List.flatten l2)) in
    Labels.add_list_of_list_pairs labels ll)

(*let rec count_multiply_list = function
    [] -> 1
  | l :: ll ->
     let x = count_multiply_list ll in
     if x > 1000000000 then x else
     Xlist.size l * x*)

let match_reaction level r labels =
(*   print_endline "match_reaction 1"; *)
  if Sys.time () -. !Types.time > Types.timeout then ((*print_endline "TIMEOUT";*) raise Timeout) else
(*   print_endline "match_reaction 2"; *)
  let hash = Hash.make_hash Hash.hash_key level r.graph labels in
(*   print_endline "match_reaction 3"; *)
  let comp = match_isomorphic_connected_components_simple r.graph r.reactant_ids r.product_ids hash level in
(*   print_endline "match_reaction 4"; *)
  if Collection.is_empty comp then Collection.empty else (
(*   print_endline "match_reaction 5"; *)
  let no_divisions =
      try Collection.int_log (Collection.number_of_product_permutations comp)
      with Collection.Empty -> raise (SolutionNotFound "no divisions") in
    if no_divisions > 3 then raise (Ambiguity "too many divisions") else
  let divisions = Collection.generate_product_permutations comp in (*!!!!*)
(*   let divisions = Collection.map comp (fun (la,lb) -> Collection.combine la lb (fun _ _ -> true)) in *)
(*   print_endline "match_reaction 6"; *)
  Collection.flatten_map divisions (fun division ->
(*     print_endline "match_reaction a"; *)
    if Sys.time () -. !Types.time > Types.timeout then ((*print_endline "TIMEOUT";*) raise Timeout) else
    let matchings = Collection.map division (fun (selection1,selection2) ->
      Isomorphism.find_isomorphisms false level hash hash r.graph r.graph selection1 selection2) in
    if Collection.int_log (Collection.number_of_products matchings) > 6 then raise (Ambiguity "too many matchings");
    Collection.flatten_product matchings))
(*  let divisions = Xlist.map (Collection.to_list divisions) Collection.to_list in
  Xlist.fold divisions [] (fun found division ->
    let matchings = Xlist.rev_map division (fun (selection1,selection2) ->
      Isomorphism.find_isomorphisms false level hash hash r.graph r.graph selection1 selection2) in
    (Xlist.rev_map (Xlist.multiply_list matchings) List.flatten) @ found)*)

let match_reaction_full_key level r labels =
(*   print_endline "match_reaction 1"; *)
  if Sys.time () -. !Types.time > Types.timeout then ((*print_endline "TIMEOUT";*) raise Timeout) else
(*   print_endline "match_reaction 2"; *)
  let hash = Hash.make_hash Hash.very_full_hash_key level r.graph labels in
(*   print_endline "match_reaction 3"; *)
  let comp = match_isomorphic_connected_components_simple r.graph r.reactant_ids r.product_ids hash level in
(*   print_endline "match_reaction 4"; *)
  if Collection.is_empty comp then Collection.empty else (
(*   print_endline "match_reaction 5"; *)
  let no_divisions =
      try Collection.int_log (Collection.number_of_product_permutations comp)
      with Collection.Empty -> raise (SolutionNotFound "no divisions") in
    if no_divisions > 3 then raise (Ambiguity "too many divisions") else
  let divisions = Collection.generate_product_permutations comp in (*!!!!*)
(*   let divisions = Collection.map comp (fun (la,lb) -> Collection.combine la lb (fun _ _ -> true)) in *)
(*   print_endline "match_reaction 6"; *)
  Collection.flatten_map divisions (fun division ->
(*     print_endline "match_reaction a"; *)
    if Sys.time () -. !Types.time > Types.timeout then ((*print_endline "TIMEOUT";*) raise Timeout) else
    let matchings = Collection.map division (fun (selection1,selection2) ->
      Isomorphism.find_isomorphisms false level hash hash r.graph r.graph selection1 selection2) in
    if Collection.int_log (Collection.number_of_products matchings) > 6 then raise (Ambiguity "too many matchings");
    Collection.flatten_product matchings))

(*let rec isomorphic_partition2 level l =
  let l = Collection.of_list l in
  let ll = Collection.quotient l (fun x y -> Isomorphism.are_isomorphic_lists true level x y) in
  Xlist.map (Collection.to_list ll) (fun l -> Xlist.map (Collection.to_list l) fst)

(*let disambiguate_matchings3 r subtasks =
  let subtasks = Xlist.map subtasks (fun (matchings,b) ->
    let reactions = Xlist.map matchings (fun (matching(*,l*)) ->
      let ll = Xlist.map matching (fun (v,w) -> [v;w]) in
      let labels = Labels.add_list_of_lists r.empty_labels ll in
      let hash = Hash.make_hash Hash.full_hash_key 1 r.graph labels in
      (matching,b),Xlist.map r.molecules (fun molecule -> hash,r.graph,molecule)) in
    let label_buckets = isomorphic_partition2 1 reactions in
    Xlist.map label_buckets List.hd) in
  List.flatten alt_labels_list*)

let rec sort_molecules_rec rev (l,molecules) (v,w) =
  if molecules = [] then l,rev else
  if IntSet.mem (List.hd molecules).ids w then (List.hd molecules) :: l, rev @ List.tl molecules
  else sort_molecules_rec ((List.hd molecules) :: rev) (l,List.tl molecules) (v,w)

let fst_compare x y = compare (fst x) (fst y)

let sort_molecules matching molecules =
  let l,molecules = Xlist.fold (List.sort fst_compare matching) ([],molecules) (sort_molecules_rec []) in
  List.rev (molecules @ l)

(*let string_of_col col =
  String.concat "; " (Xlist.map (Collection.to_list col) (fun (x,y) ->
    Printf.sprintf "%d,%d" x y))*)

let disambiguate_matchings3 r subtasks =
  let subtasks = Xlist.map subtasks (fun (cand,matching,hydrogens) ->
(*          Printf.printf "disambiguate_matchings3 \n  cand=%s \n  matching=%s \n  hydrogens=%s\n%!" (string_of_col cand) (string_of_col matching)
           (String.concat ";" (Xlist.map (IntSet.to_list hydrogens) string_of_int));*)
      let ll = Xlist.fold (Collection.to_list matching) [IntSet.to_list hydrogens] (fun ll (v,w) -> if IntSet.mem hydrogens v then ll else [v;w] :: ll) in
      let labels = Labels.add_list_of_lists r.empty_labels ll in
      let hash = Hash.make_hash Hash.full_hash_key 1 r.graph labels in
      (cand,matching,hydrogens),Xlist.map (r.reactants @ (sort_molecules (Collection.to_list matching) r.products)) (fun m -> hash,r.graph,Collection.of_list (IntSet.to_list m.ids))) in
  let label_buckets = isomorphic_partition2 1 subtasks in
  Xlist.map label_buckets List.hd*)

let disambiguate_matchings r subtasks =
  let subtasks = Collection.of_list subtasks in
  let subtasks_ann = Collection.annotate subtasks (fun (cand,matching,hydrogens) ->
    let ll = Xlist.fold (Collection.to_list matching) [IntSet.to_list hydrogens] (fun ll (v,w) -> if IntSet.mem hydrogens v then ll else [v;w] :: ll) in
    let labels = Labels.add_list_of_lists r.empty_labels ll in
    r.graph, Hash.make_hash Hash.full_hash_key 1 r.graph labels) in
  let label_buckets =
    try Collection.quotient subtasks_ann (fun (x1,(graph1,hash1)) (x2,(graph2,hash2)) ->
      are_isomorphic_graphs graph1 graph2 r.reactant_ids r.reactant_ids hash1 hash2 1 &&
      are_isomorphic_graphs graph1 graph2 r.product_ids r.product_ids hash1 hash2 1)
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during partitioning") in
  Collection.to_list (Collection.deannotate (Collection.reverse_projection label_buckets))

(*let disambiguate_matchings r subtasks =
  let subtasks = Collection.of_list subtasks in
  let subtasks_ann = Collection.annotate subtasks (fun (cand,matching,hydrogens) ->
    let ll = Xlist.fold (Collection.to_list matching) [IntSet.to_list hydrogens] (fun ll (v,w) -> if IntSet.mem hydrogens v then ll else [v;w] :: ll) in
    let labels = Labels.add_list_of_lists r.empty_labels ll in
    r.graph, Hash.make_hash Hash.full_hash_key 1 r.graph labels) in
  let label_buckets =
    try Collection.quotient subtasks_ann (fun x y -> Isomorphism.are_isomorphic_alt_lists true 1 x y (Collection.to_list r.component_perms))
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during partitioning") in
  Collection.to_list (Collection.deannotate (Collection.reverse_projection label_buckets))*)

(*let print_graph g h =
  Int.iter 0 (Array.length g-1) (fun i ->
    let p,l = g.(i) in
    if p.name <> "H" then printf "%s_%d %s: %s\n%!" p.name p.id h.(i) (String.concat " " (Xlist.map l (fun (b,p) -> Hash.key_of_bond b ^ string_of_int p.id))))*)

(*let disambiguate_aroma subtasks =
  if subtasks = [] then [] else
  let component_perms = try Collection.to_list (List.hd subtasks).component_perms with _ -> failwith "disambiguate_aroma" in
  let subtasks = Collection.of_list subtasks in
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

let disambiguate_rearrangements subtasks =
  if subtasks = [] then [] else
  let subtasks = Collection.of_list subtasks in
  let subtasks_ann = Collection.annotate subtasks (fun r ->
    r.graph, r.reactant_ids, r.product_ids, Hash.make_hash Hash.full_hash_key 1 r.graph r.empty_labels) in
  let label_buckets =
    try Collection.quotient subtasks_ann (fun (x1,(graph1,reactant_ids1,product_ids1,hash1)) (x2,(graph2,reactant_ids2,product_ids2,hash2)) ->
      are_isomorphic_graphs graph1 graph2 reactant_ids1 reactant_ids2 hash1 hash2 1 &&
      are_isomorphic_graphs graph1 graph2 product_ids1 product_ids2 hash1 hash2 1)
    with Isomorphism.TreeToBig -> raise (Ambiguity "decision tree to big during partitioning") in
  Collection.to_list (Collection.deannotate (Collection.reverse_projection label_buckets))

(*let disambiguate_rearrangements subtasks =(*subtasks*)
  if subtasks = [] then [] else
  let subtasks = Collection.annotate (Collection.of_list subtasks) (fun r ->
    let hash = (Hash.make_hash Hash.full_hash_key 3 r.graph r.empty_labels).(3) in
    let conn_r = get_connected_components Collection.empty r.graph r.reactant_ids in
    let conn_p = get_connected_components Collection.empty r.graph r.product_ids in
    let conn_r = Collection.map conn_r (fun c -> String.concat " " (Collection.to_list (Collection.sort (Collection.map c (fun v -> hash.(v)))))) in
    let conn_p = Collection.map conn_p (fun c -> String.concat " " (Collection.to_list (Collection.sort (Collection.map c (fun v -> hash.(v)))))) in
    String.concat "." (Collection.to_list conn_r) ^ ">>" ^ String.concat "." (Collection.to_list conn_p)) in
  let subtasks = Collection.group2 subtasks in
  Collection.to_list (Collection.reverse_projection (Collection.deannotate subtasks))*)

(*********************************************************************************************)

type history = A of int | N of int | AB of int | NB of int | P of int | B of int

let string_of_history history =
  String.concat "  " (Xlist.map history (function
      A n -> "A " ^ string_of_int n
    | AB n -> "AB " ^ string_of_int n
    | B n -> "B " ^ string_of_int n
    | N n -> "N " ^ string_of_int n
    | NB n -> "NB " ^ string_of_int n
    | P n -> "P " ^ string_of_int n))

type action = Parse of int | BorderParse of int | Parsed of int | ParsedAmb of int | FailAmb | FailNone | Unknown

let get_level3 = function
    [] -> Parse 4
  | P 1 :: _ -> Parsed 1
  | P n :: _ -> BorderParse (n-1)
  | B 1 :: _ -> Parsed 1
  | B n :: _ -> Parse n
  | N 4 :: _ -> Parse 3
  | N 1 :: _ -> FailNone
  | N n :: _ -> BorderParse (n-1)
  | A _ :: _ -> FailAmb
  | NB n :: _ -> Parse n
  | AB _ :: _ -> FailAmb
(*   | _ -> Unknown *)

let rec multilevel_label_reaction3 r history labels_list =
  match get_level3 history with
    Parse level ->
      (try
        let alt_labels_list =
          if labels_list = [] then disambiguate_labels level r (label_reaction level r r.empty_labels)
          else List.flatten (Xlist.map labels_list (fun labels ->
            try disambiguate_labels level r (label_reaction level r labels) with SolutionNotFound _ -> [])) in
        if alt_labels_list = [] then raise (SolutionNotFound "cummulative") else
        multilevel_label_reaction3 r ((P level) :: history) alt_labels_list
      with Ambiguity _ -> multilevel_label_reaction3 r ((A level) :: history) labels_list
         | SolutionNotFound _ -> multilevel_label_reaction3 r ((N level) :: history) labels_list)
  | BorderParse level ->
      (try
        let alt_labels_list =
          if labels_list = [] then disambiguate_labels level r (label_border_reaction level r r.empty_labels)
          else List.flatten (Xlist.map labels_list (fun labels ->
            try disambiguate_labels level r (label_border_reaction level r labels) with SolutionNotFound _ -> [])) in
        if alt_labels_list = [] then raise (SolutionNotFound "cummulative") else
        multilevel_label_reaction3 r ((B level) :: history) alt_labels_list
      with Ambiguity _ -> multilevel_label_reaction3 r ((AB level) :: history) labels_list
         | SolutionNotFound _ -> multilevel_label_reaction3 r ((NB level) :: history) labels_list)
  | result -> result, List.rev history, labels_list


(**
let rec has_pattern_rec ma mb = function
    [] -> false
  | x :: l ->
      (match create_decision_tree_pat ma mb IntMap.empty IntMap.empty x with
        Tfinish -> true
      | Terror -> has_pattern_rec ma mb l
      | _ -> failwith "has_pattern_rec")

let string_of_graph m =
  String.concat "\n" (IntMap.fold m [] (fun t n (s,l) ->
    (Printf.sprintf "%s %d: %s" s n (String.concat " " (Xlist.map l (fun (s,n) -> s ^ " " ^ string_of_int n)))) :: t))

let has_pattern a b (*smile*) =
  match (remove_root_ions a,remove_root_ions b) with
    SAtom(sa,pa,[],na),SAtom(sb,pb,[],nb) -> sa = sb
  | _ -> (
  let ca = SummaricModel.count_atoms StringQMap.empty a in
(*   Printf.printf "SAtom count: %s\n%!" (string_of_counted_atoms ca); *)
  let cb = SummaricModel.count_atoms StringQMap.empty b in
  let c = SummaricModel.compare_atom_counts ca cb in
  if not (SummaricModel.is_submolecule c) then false else
(*   Printf.printf "Compared atom counts: [%s]\n%!" (String.concat ";" (Xlist.map c (fun (s,v1,v2) -> Printf.sprintf "%s,%d,%d" s v1 v2))); *)
  let s = select_matching_root c in
(*   Printf.printf "Matching root: '%s'\n%!" s; *)
  let ma = make_atom_graph a in
  let mb = make_atom_graph b in
  let l = (*try*) create_starting_assumptions_pat s ma mb in
    (*with Not_found ->*)
(*      Printf.printf "Pattern graph: %s\n%!" (string_of_graph ma);
      Printf.printf "Molecule smile: %s\n%!" smile;
      Printf.printf "Molecule graph: %s\n%!" (string_of_graph mb);
      Printf.printf "SAtom count a: %s\n%!" (SummaricModel.string_of_counted_atoms ca);
      Printf.printf "SAtom count b: %s\n%!" (SummaricModel.string_of_counted_atoms cb);
      Printf.printf "Compared atom counts: [%s]\n%!" (String.concat ";" (Xlist.map c (fun (s,v1,v2) -> Printf.sprintf "%s,%d,%d" s v1 v2)));
      Printf.printf "Matching root: '%s'\n%!" s; *)
    (*  failwith "has_pattern"
    in  *)
  counter := 0;
  has_pattern_rec ma mb l)

let test_has_pattern patterns filename =
  let time1 = Unix.gettimeofday () in
(*   let molecules = Import.load_molecules () in *)
  let molecules = Import.load_molecules () in
  let time2 = Unix.gettimeofday () in
  Printf.printf "Loading time: %f\n" (time2 -. time1);
  Printf.printf "No molecules: %d\n" (Xlist.size molecules);
  Xlist.iter patterns (fun pat ->
    Printf.printf "\nLooking for %s\n%!" pat;
    let time2 = Unix.gettimeofday () in
    let pat = Smiles.parse_smile_molecule "x" pat in
    let pat,next_id = Smiles.assign_unique_ids 1 pat in
    let solution = List.rev (Xlist.rev_map molecules (fun s ->
(*     print_endline s; *)
      let t = Smiles.parse_smile_molecule "x" s in
      let t,_ = Smiles.assign_unique_ids next_id t in
      let b = has_pattern pat t (*s*) in
(*       if b then print_endline ("found: " ^ s); *)
      b)) in
    let time3 = Unix.gettimeofday () in
    File.file_out filename (fun file ->
      Xlist.iter solution (fun b ->
        if b then Printf.fprintf file "Y\n" else Printf.fprintf file "N\n"));
    let time4 = Unix.gettimeofday () in
    Printf.printf "Matching time: %f\n" (time3 -. time2);
    Printf.printf "Avg matching time of a molecule: %f\n" ((time3 -. time2) /. (float (Xlist.size molecules)));
    Printf.printf "No found molecules: %d\n" (Xlist.fold solution 0 (fun n b -> if b then n+1 else n));
    Printf.printf "Writing time: %f\n" (time4 -. time3);
    () )

let patterns = [
  "O=C1CCN1";
  "CC(C(O)=O)c1ccccc1";
  "OC(=O)C1=CSC2CC(=O)N12";
  "C1CC2CCC3C(CCC4CCCCC34)C2C1";
  "CN1CCC23C4CCCC2C1Cc1ccc(O)c(O4)c31";
  "c1ccccc1";
  ]

(** testowanie ile czasu trwa wyszukanie wzorca *)
(* let _ =  test_has_pattern patterns "results/matching.txt" *)
**)
