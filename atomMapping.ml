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

(*let string_of_col col =
  String.concat "; " (Xlist.map (Collection.to_list col) (fun (x,y) ->
    Printf.sprintf "%d,%d" x y))*)

let find_candidate_bonds graph ids labels (*unbreakable*) =
  Collection.of_list (Int.fold 1 (Array.length graph - 1) [] (fun bonds i ->
    if IntSet.mem ids i then
      Xlist.fold (snd graph.(i)) bonds (fun bonds (_,q) ->
        if (Labels.mem labels i && Labels.mem labels q.id) || q.id < i (*|| Xlist.mem unbreakable (i,q.id)*) then bonds else (i,q.id) :: bonds)
    else bonds))

let create_candidates r labels n =
(*   printf "create_candidates 1\n%!";  *)
  let reactant_bonds = find_candidate_bonds r.graph r.reactant_ids labels (*r.unbreakable*) in
  if Collection.size reactant_bonds < n then Collection.empty else (
(*   printf "create_candidates 2\n%!"; *)
  let product_bonds = find_candidate_bonds r.graph r.product_ids labels (*[]*) in
(*   printf "create_candidates 3\n%!"; *)
  let x = Collection.group (Collection.flatten_map (Collection.generate_partition_combinations n reactant_bonds) (fun (ex,reactant_bonds) ->
(*     printf "create_candidates 3a\n%!";  *)
    let bond_buckets = Collection.string_pair_bucket_group reactant_bonds product_bonds (fun (i,j) ->
      String.concat "-" (List.sort compare [(fst r.graph.(i)).name; (fst r.graph.(j)).name])) in
(*     printf "create_candidates 3b %d\n%!" (Collection.int_log (Collection.number_of_exclusion_schemes bond_buckets)); *)
    let no_ex_sch =
      try Collection.int_log (Collection.number_of_exclusion_schemes bond_buckets)
      with Collection.Empty -> 0 in
    if no_ex_sch > 6 then raise (CommonSubstructure.Ambiguity "too many exclusion schemes for candidates") else
    let excl = Collection.get_exclusion_schemas bond_buckets in
(*     printf "create_candidates 3c\n%!"; *)
    let excl = Collection.map excl (fun e -> Collection.sort (Collection.sum ex e), labels) in
(*     printf "create_candidates 3d\n%!"; *)
    excl)) in
(*   printf "create_candidates 4\n%!"; *)
(*  Xlist.iter (Collection.to_list x) (fun (e,_) ->
    printf "e=%s\n%!" (String.concat " " (Xlist.map (Collection.to_list e) (fun (i,j) -> sprintf "%d-%d" i j))));*)
  let x = Collection.flatten_map x (fun (e,labels) ->
(*     printf "create_candidates 4a\n%!"; *)
    if Collection.is_empty (Collection.intersection e r.unbreakable) then (
(*       printf "e3=%s\n%!" (String.concat " " (Xlist.map (Collection.to_list e) (fun (i,j) -> sprintf "%d-%d" i j))); *)
(*       printf "create_candidates 4b\n%!"; *)
      Collection.singleton (e,labels) ) else (
(*       printf "e4=%s\n%!" (String.concat " " (Xlist.map (Collection.to_list e) (fun (i,j) -> sprintf "%d-%d" i j))); *)
(*       printf "e5=%s\n%!" (String.concat " " (Xlist.map (Collection.to_list (Collection.intersection e r.unbreakable)) (fun (i,j) -> sprintf "%d-%d" i j))); *)
(*       printf "create_candidates 4c\n%!";  *)
      Collection.empty)) in
(*  Xlist.iter (Collection.to_list x) (fun (e,_) ->
    printf "e2=%s\n%!" (String.concat " " (Xlist.map (Collection.to_list e) (fun (i,j) -> sprintf "%d-%d" i j))));*)
(*   printf "create_candidates 5\n%!";  *)
  x)

let max_global_broken_bonds = 1000

let stratify_candidates cands =
  let a = Array.make max_global_broken_bonds [] in
  let size = Xlist.fold (Collection.to_list cands) 0 (fun size (bonds,labels) ->
    let i = Collection.size bonds in
    if i >= max_global_broken_bonds then failwith "stratify_candidates" else
    a.(i) <- (bonds,labels) :: a.(i);
    max i size) in
  a,size

let remove_subset_cands good_cands cands = (* FIXME: nie ma scalania identycznych kandydatów *)
  Xlist.fold (Collection.to_list good_cands) cands (fun cands (good_cand,_,_) -> (* FIXME: czy na pewno taka kolejność argumentów ? *)
    Collection.of_list (Xlist.fold (Collection.to_list cands) [] (fun cands (cand,labels) ->
      if Collection.size good_cand = Collection.size (Collection.intersection good_cand cand) then cands else (cand,labels) :: cands)))

let remove_bond2a g v w =
  let p,l = g.(v) in
  let l = Xlist.fold l [] (fun l (b,q) -> if q.id = w then l else (b,q) :: l) in
  g.(v) <- p,l

let rec remove_list w = function
    [] -> []
  | v :: l -> if v = w then (*remove_list w*) l else v :: remove_list w l

let remove_bond2a_pat g v w =
  let p,l = g.(v) in
  let l = Xlist.fold l [] (fun l (b,q) -> if q.id = w then l else (b,q) :: l) in
  g.(v) <- {p with hydrogens=remove_list w p.hydrogens; fluors=remove_list w p.fluors},l

let remove_cand graph cand =
  let g = Array.copy graph in
  Xlist.iter (Collection.to_list cand) (fun (v,w) ->
        remove_bond2a g v w;
        remove_bond2a g w v);
  g

let remove_cand_pat graph cand =
  let g = Array.copy graph in
  Xlist.iter (Collection.to_list cand) (fun (v,w) ->
        remove_bond2a_pat g v w;
        remove_bond2a_pat g w v);
  g

let rec match_hydrogens_and_fluors_rec (matching,matched_hf_reactants,matched_hf_products) = function (* FIXME: może zmaczować H z F *)
    [],_ -> matching,matched_hf_reactants,matched_hf_products
  | _,[] -> matching,matched_hf_reactants,matched_hf_products
  | i :: l1, j :: l2 -> match_hydrogens_and_fluors_rec ((i,j) :: matching, IntSet.add matched_hf_reactants i, IntSet.add matched_hf_products j) (l1,l2)

(* FIXME: usunięcie wiązań z wodorami wskazywanymi przez szablon powoduje, że wodory te są traktowane jako odrywane nawet, gdy faktycznie nie ma to miejsca *)
let match_hydrogens_and_fluors cand r matching =
  let matching,matched_hf_reactants,matched_hf_products = Xlist.fold matching (matching,IntSet.empty,IntSet.empty) (fun found (i,j) ->
    let found = match_hydrogens_and_fluors_rec found ((fst r.graph.(i)).hydrogens,(fst r.graph.(j)).hydrogens) in
    match_hydrogens_and_fluors_rec found ((fst r.graph.(i)).fluors,(fst r.graph.(j)).fluors)) in
  let lf1 = Xlist.fold r.reactants [] (fun l m -> IntSet.fold m.hf_ids l (fun l id -> if IntSet.mem matched_hf_reactants id || (fst r.graph.(id)).name <> "F" then l else id :: l)) in
  let lf2 = Xlist.fold r.products [] (fun l m -> IntSet.fold m.hf_ids l (fun l id -> if IntSet.mem matched_hf_products id || (fst r.graph.(id)).name <> "F" then l else id :: l)) in
  let llf = Pair.make_matchings (Pair.singleton lf1 lf2) in
  let matchings = Xlist.map llf (fun l -> l @ matching) in
  let lh1 = Xlist.fold r.reactants [] (fun l m -> IntSet.fold m.hf_ids l (fun l id -> if IntSet.mem matched_hf_reactants id || (fst r.graph.(id)).name <> "H" then l else id :: l)) in
  let lh2 = Xlist.fold r.products [] (fun l m -> IntSet.fold m.hf_ids l (fun l id -> if IntSet.mem matched_hf_products id || (fst r.graph.(id)).name <> "H" then l else id :: l)) in
(*  let llh = Pair.make_matchings (Pair.singleton lh1 lh2) in
  List.flatten (Xlist.rev_map matchings (fun matching -> Xlist.map llh (fun l -> l @ matching)))*)
  (* print_endline (String.concat " " (Xlist.map lh1 string_of_int));
  print_endline (String.concat " " (Xlist.map lh2 string_of_int)); *)
  let lhh = try List.combine lh1 lh2 with _ -> failwith "match_hydrogens_and_fluors: invalid stoichimetry" in
(*   let lhh = Xlist.map lh1 (fun h -> h,0) @ Xlist.map lh2 (fun h -> 0,h) in *)
  let hydrogens = IntSet.of_list (lh1 @ lh2) in
  Collection.of_list (Xlist.rev_map matchings (fun matching -> cand,Collection.of_list (lhh @ matching),hydrogens))

let select_carbonyls r matchings =
  Collection.of_list (List.rev (Xlist.fold (Collection.to_list matchings) [] (fun matchings matching ->
    let b = Xlist.fold matching true (fun b (i,j) ->
      match r.graph.(i),r.graph.(j) with
        ({name="O"},[Double,_]),(_,[]) -> b
      | ({name="O"},[Double,{id=c}]),(_,l) ->
           if (fst r.graph.(c)).name <> "C" then b else
           let d = try Xlist.assoc matching c with Not_found -> failwith "select_carbonyls" in
           let f = Xlist.fold l false (fun f (_,p) -> if p.id = d then true else f) in
           if f then b else ((*Printf.printf "%d-%d\n" i j;*) false)
      | _ -> b) in
    if b then matching :: matchings else matchings)))

let rec create_good_candidates_rec msg r cands_array i cands_array_size (good_cands:((int * int) Collection.collection * (int * int) Collection.collection * IntSet.t) Collection.collection) =
  if i > cands_array_size then good_cands else
  if Sys.time () -. !Types.time > Types.timeout then ((*print_endline "TIMEOUT";*) raise Timeout) else
  let cands = remove_subset_cands good_cands (Collection.of_list cands_array.(i)) in
(*   let time1 = Sys.time () in
  printf "create_good_candidates_rec 1 |cands|=%d\n%!" (Collection.size cands); *)
  let cands:((int * int) Collection.collection * (int * int) Collection.collection * IntSet.t) Collection.collection = Collection.flatten_map cands (fun (cand,l) ->
    let r2 =  {r with graph = remove_cand r.graph cand} in
    Collection.flatten_map l (fun labels ->
      try
(*         print_endline "create_good_candidates_rec a"; *)
        let matchings = CommonSubstructure.match_reaction 4 r2 labels in
(*         let matchings = select_carbonyls r matchings in  (* !!! wyłączone testowo do OrgSyn *)  *)
(*         print_endline "create_good_candidates_rec b"; *)
        let x = Collection.flatten_map matchings (fun matching ->
          match_hydrogens_and_fluors cand r2 (Collection.to_list matching)) in
(*         print_endline "create_good_candidates_rec c"; *)
        x
(*
        Collection.of_list in
(*                Xlist.fold matchings good_cands (fun good_cands matching ->
                  let matchings, hydrogens = match_hydrogens_and_fluors r matching in
                  Xlist.fold matchings good_cands (fun good_cands matching -> (cand,labels,bonds,matching,hydrogens) :: good_cands))*)
        Collection.annotate matchings (fun matching -> cand)*)
      with Timeout -> raise Timeout
      | e -> (msg := !msg @ [Printexc.to_string e]; Collection.empty))) in
(*   let time2 = Sys.time () in
  printf "create_good_candidates_rec 2 time=%.4f\n%!" (time2 -. time1);  *)
  let good_cands = Collection.sum cands good_cands in
  create_good_candidates_rec msg r cands_array (i+1) cands_array_size good_cands

(*let rec create_good_candidates_rec msg r cands_array i cands_array_size (good_cands:((int * int) Collection.collection * (int * int) Collection.collection * IntSet.t) Collection.collection) =
  if i > cands_array_size then good_cands else
  let cands = remove_subset_cands good_cands (Collection.of_list cands_array.(i)) in
  let time1 = Unix.gettimeofday () in
  let cands:((int * int) Collection.collection * (int * int) Collection.collection * IntSet.t) Collection.collection =
    let work = Collection.map cands (fun (cand,l) ->
      let r2 =  {r with graph = remove_cand r.graph cand} in
      r,cand,r2,l) in
    let s = Marshal.to_string work [Marshal.No_sharing] in
    printf "create_good_candidates_rec 1 |cands|=%d |work|=%d\n%!" (Collection.size cands) (String.length s);
    let _,work = Xlist.fold (Collection.to_list work) (1,[]) (fun (id,work) t ->
      id+1, (string_of_int id, t) :: work) in
    Collection.flatten (Collection.of_list (Overseer.execution Overseer.io_list msg work)) in
  let time2 = Unix.gettimeofday () in
  printf "create_good_candidates_rec 2 time=%.4f\n%!" (time2 -. time1);
  let good_cands = Collection.sum cands good_cands in
  create_good_candidates_rec msg r cands_array (i+1) cands_array_size good_cands*)

let rec create_good_candidates msg max_broken_bonds no_broken_bonds r alt_labels_list =
(*   print_endline "create_good_candidates 1"; *)
  let cands = Collection.flatten_map alt_labels_list (fun labels ->
    create_candidates r labels no_broken_bonds) in
(*   print_endline "create_good_candidates 2"; *)
  if Collection.is_empty cands then no_broken_bonds,cands,Collection.empty else (
(*   print_endline "create_good_candidates 3"; *)
  msg := !msg @ [(*Printf.sprintf "no_broken_bonds=%d no candidates=%d" no_broken_bonds (Collection.size cands)*)];
  if Collection.size cands > 1000000 then (
    msg := !msg @ ["TO MANY CANDIDATES"];
    no_broken_bonds,cands,Collection.empty) else (
(*   print_endline "create_good_candidates 4"; *)
  let cands_array,cands_array_size = stratify_candidates cands in
  let (good_cands:'a Collection.collection) = create_good_candidates_rec msg r cands_array 0 cands_array_size Collection.empty in
(*   print_endline "create_good_candidates 5"; *)
  if Collection.is_empty good_cands && no_broken_bonds < max_broken_bonds then create_good_candidates msg max_broken_bonds (no_broken_bonds+1) r alt_labels_list
  else no_broken_bonds,cands,(good_cands:'a (*(int * int) Collection.collection * (int * int) Collection.collection*) Collection.collection)))

let rec create_good_candidates2 msg max_broken_bonds no_broken_bonds arl =
(*   print_endline "create_good_candidates2 1"; *)
(*   printf "no_broken_bonds=%d\n%!" no_broken_bonds;  *)
  let rcands = Xlist.fold arl [] (fun rcands (r,alt_labels_list) ->
(*     printf "create_good_candidates2 a %d\n%!" (Collection.size alt_labels_list); *)
(*     printf "unbreakable=%s\n%!" (String.concat " " (Xlist.map (Collection.to_list r.unbreakable) (fun (i,j) -> sprintf "%d-%d" i j))); *)
    let cands = Collection.flatten_map alt_labels_list (fun labels ->
(*       print_endline "create_good_candidates2 c"; *)
      create_candidates r labels no_broken_bonds) in
(*     print_endline "create_good_candidates2 b"; *)
    if Collection.is_empty cands then ((*print_endline "empty cands";*) rcands) else (r,cands) :: rcands) in
(*   print_endline "create_good_candidates2 2"; *)
(*   printf "|cands|=%d\n%!" (Xlist.fold rcands 0 (fun n (_,cands) -> n + Collection.size cands)); *)
  msg := !msg @ [(*Printf.sprintf "no_broken_bonds=%d no candidates=%d" no_broken_bonds (Collection.size cands)*)];
  if rcands = [] then no_broken_bonds,Collection.empty,[] else
  if Xlist.fold rcands 0 (fun n (_,cands) -> n + Collection.size cands) > 1000000 then (
    msg := !msg @ ["TO MANY CANDIDATES"];
    no_broken_bonds,Collection.empty(*cands*),[]) else (
(*   print_endline "create_good_candidates2 3"; *)
  let rgood_cands = Xlist.fold rcands [] (fun rgood_cands (r,cands) ->
    let cands_array,cands_array_size = stratify_candidates cands in
    let (good_cands:'a Collection.collection) = create_good_candidates_rec msg r cands_array 0 cands_array_size Collection.empty in
    if Collection.is_empty good_cands then rgood_cands else (r,good_cands) :: rgood_cands) in
(*   print_endline "create_good_candidates2 4"; *)
  if rgood_cands = [] && no_broken_bonds < max_broken_bonds then create_good_candidates2 msg max_broken_bonds (no_broken_bonds+1) arl
  else no_broken_bonds,Collection.empty(*cands*),rgood_cands)

(*****************************)

let rec set_invisible_hydrogens visible rev = function
    Atom(p,l) ->
       let rev = if p.name = "H" && not (IntSet.mem visible p.id) then p.id :: rev else rev in
       Xlist.fold l rev (fun rev (_,a) -> set_invisible_hydrogens visible rev a)
  | Link _ -> rev

let find_reaction_center_bonds r matching =
  let matching = Xlist.fold matching IntMap.empty (fun map (v,w) ->
    let map = IntMap.add map v w in
    IntMap.add map w v) in
  let center1,center2 = Int.fold 1 (Array.length r.graph - 1) ([],[]) (fun (center1,center2) i ->
    if not (IntMap.mem matching i) then center1,center2 else
    let l1 = Xlist.fold (snd r.graph.(i)) [] (fun l1 (b,q) -> if IntMap.mem matching q.id then (b, {q with id = IntMap.find matching q.id}) :: l1 else l1) in
    let l2 = snd r.graph.(IntMap.find matching i) in
    let ll = Pair.string_bucket_group (Pair.singleton l1 l2) (fun (b,q) ->
              (* if q.name = "H" then raise Not_found else*) Hash.key_of_bond b ^ string_of_int (IntMap.find matching q.id)) (fun (_,q) -> q.id) in
    let center2 = Pair.fold ll center2 (fun center -> function
        [_],[_] -> center
      | [],[id] -> (IntMap.find matching i, id) :: center
      | [id],[] -> (i, IntMap.find matching id) :: center
      | _ -> failwith "find_reaction_center_bonds") in
    let ll = Pair.string_bucket_group (Pair.singleton l1 l2) (fun (b,q) ->
               (*if q.name = "H" then raise Not_found else*) string_of_int (IntMap.find matching q.id)) (fun (_,q) -> q.id) in
    let center1 = Pair.fold ll center1 (fun center -> function
        [_],[_] -> center
      | [],[id] -> (IntMap.find matching i, id) :: center
      | [id],[] -> (i, IntMap.find matching id) :: center
      | _ -> failwith "find_reaction_center_bonds") in
    center1,center2) in
  let center_atoms = Xlist.fold center2 IntSet.empty (fun center_atoms (i,j) ->
    IntSet.add (IntSet.add center_atoms i) j) in
  let center = Xlist.fold center2 IntMap.empty (fun center (i,j) ->
      let v = min i j in
      let w = max i j in
      IntMap.add_inc center v (IntMap.add IntMap.empty w 2) (fun map -> IntMap.add map w 2)) in
  let center = Xlist.fold center1 center (fun center (i,j) ->
      let v = min i j in
      let w = max i j in
      IntMap.add_inc center v (IntMap.add IntMap.empty w 1) (fun map -> IntMap.add map w 1)) in
  center, center_atoms

let find_reaction_center_bonds2 r matching =
  let matching = Xlist.fold matching IntMap.empty (fun map (v,w) ->
    let map = IntMap.add map v w in
    IntMap.add map w v) in
  let center1,center2 = Int.fold 1 (Array.length r.original_graph - 1) ([],[]) (fun (center1,center2) i ->
    if not (IntMap.mem matching i) then center1,center2 else
    let l1 = Xlist.fold (snd r.original_graph.(i)) [] (fun l1 (b,q) -> if IntMap.mem matching q.id then (b, {q with id = IntMap.find matching q.id}) :: l1 else l1) in
    let l2 = snd r.original_graph.(IntMap.find matching i) in
    let ll = Pair.string_bucket_group (Pair.singleton l1 l2) (fun (b,q) ->
              (* if q.name = "H" then raise Not_found else*) Hash.key_of_bond b ^ string_of_int (IntMap.find matching q.id)) (fun (_,q) -> q.id) in
    let center2 = Pair.fold ll center2 (fun center -> function
        [_],[_] -> center
      | [],[id] -> (IntMap.find matching i, id) :: center
      | [id],[] -> (i, IntMap.find matching id) :: center
      | _ -> failwith "find_reaction_center_bonds") in
    let ll = Pair.string_bucket_group (Pair.singleton l1 l2) (fun (b,q) ->
               (*if q.name = "H" then raise Not_found else*) string_of_int (IntMap.find matching q.id)) (fun (_,q) -> q.id) in
    let center1 = Pair.fold ll center1 (fun center -> function
        [_],[_] -> center
      | [],[id] -> (IntMap.find matching i, id) :: center
      | [id],[] -> (i, IntMap.find matching id) :: center
      | _ -> failwith "find_reaction_center_bonds") in
    center1,center2) in
  let center_atoms = Xlist.fold center2 IntSet.empty (fun center_atoms (i,j) ->
    IntSet.add (IntSet.add center_atoms i) j) in
  let center = Xlist.fold center2 IntMap.empty (fun center (i,j) ->
      let v = min i j in
      let w = max i j in
      IntMap.add_inc center v (IntMap.add IntMap.empty w 2) (fun map -> IntMap.add map w 2)) in
  let center = Xlist.fold center1 center (fun center (i,j) ->
      let v = min i j in
      let w = max i j in
      IntMap.add_inc center v (IntMap.add IntMap.empty w 1) (fun map -> IntMap.add map w 1)) in
  center, center_atoms
(*
module OrderedBond = struct
  type t = int * int
  let compare = compare
end

module BondSet = Xset.Make(OrderedBond)

module OrderedCandidates = struct
  type t = BondSet.t
  let compare x y = compare (List.sort compare (BondSet.to_list x)) (List.sort compare (BondSet.to_list y))
end

module CandidateSet = Xset.Make(OrderedCandidates)
module CandidateMap = Xmap.Make(OrderedCandidates)

let extend_candidate cands cand bonds =
  BondSet.fold bonds cands (fun cands bond -> CandidateSet.add cands (BondSet.add cand bond))

let extend_candidates cands bonds =
  CandidateSet.fold cands CandidateSet.empty (fun cands cand -> extend_candidate cands cand bonds)

*)


let make_reid graph matching =
  let map = Xlist.fold matching IntMap.empty (fun map (i,j) ->
    let map = IntMap.add map i i in
    IntMap.add map j i) in
  Int.fold 1 (Array.length graph - 1) map (fun map i ->
    if IntMap.mem map i then map else
    if (fst graph.(i)).name = "H" || (fst graph.(i)).name = "F" then IntMap.add map i 0 else map)


let expand_hydrogens_and_fluors graph =
  let g = Array.copy graph in
  Array.map (fun (p,l) ->
    {p with hydrogens=[]; fluors=[]},
    Xlist.map (p.hydrogens @ p.fluors) (fun id -> Single, fst graph.(id)) @ l) g

let simplify_reid labels reid =
  let set = IntMap.fold reid IntSet.empty (fun set k v ->
    if Labels.get labels k = -1 then set else IntSet.add set v) in
  let map,_ = Xlist.fold (List.sort compare (IntSet.to_list set)) (IntMap.empty,1) (fun (map,n) v ->
    IntMap.add map v n, n+1) in
  IntMap.mapi reid (fun k v ->
    if Labels.get labels k = -1 then 0 else IntMap.find map v)

let rec get_bond_value_rec j = function
    [] -> failwith "get_bond_value_rec"
  | (b,q) :: l -> if q.id = j then b else get_bond_value_rec j l

let get_bond_value r i j =
  let b = get_bond_value_rec j (snd r.graph.(i)) in
  Smiles.int_of_bond b

let select_minimal_candidates r cands =
  (*fst*) (Xlist.fold (Collection.to_list cands) ([],max_int) (fun (selected,quality) (cand,matching,(*,labels,bonds,*)hydrogens) ->
(*    Printf.printf "select_minimal_candidates \n  cand=%s \n  matching=%s \n  hydrogens=%s\n%!" (string_of_col cand) (string_of_col matching)
      (String.concat ";" (Xlist.map (IntSet.to_list hydrogens) string_of_int));*)
    let center_bonds,_ = find_reaction_center_bonds r (Collection.to_list matching) in
    let q = IntMap.fold center_bonds 0 (fun q i map ->
      if IntSet.mem r.reactant_ids i then
        IntMap.fold map q (fun q j -> function
            1 -> q + get_bond_value r i j (*q + 10000*)
          | 2 -> q + max 0 ((get_bond_value r i j) - (get_bond_value r (Xlist.assoc (Collection.to_list matching) i) (Xlist.assoc (Collection.to_list matching) j)))(*q + 1*)
          | _ -> failwith "select_minimal_candidates")
      else q) in
    if q = quality then (cand,matching,hydrogens) :: selected, quality else
    if q < quality then [cand,matching,hydrogens], q else selected,quality))

let select_minimal_candidates2 rcands =
  (*fst*) (Xlist.fold rcands ([],max_int) (fun (selected,quality) (r,cands) ->
    let cands,q = select_minimal_candidates r cands in
    if q = quality then ({r with msg=sprintf "quality=%d" q :: r.msg},cands) :: selected, quality else
    if q < quality then [{r with msg=sprintf "quality=%d" q :: r.msg},cands], q else selected,quality))

let sort_uniq_matchings r cands =
(*   Xlist.rev_map cands (fun  *)
  cands

let remove_stoi (r,labels,center_bonds_modif,center_bonds,reid) =
  (*{r with }*)r,
  Labels.add_invisible_list labels ((IntSet.to_list r.reactant_stoi) @ (IntSet.to_list r.product_stoi)),
  center_bonds_modif,
  center_bonds,
  reid

let make_msg pref l =
  if l = [] then [] else
  [pref ^ ": " ^ String.concat " " (List.rev (StringSet.to_list (StringSet.of_list l)))]

let make_id r reid i =
  if (fst r.graph.(i)).name = "H" then "H" else
  string_of_int (try IntMap.find reid i with Not_found -> i)

let center_bonds_msg r center_bonds reid =
  let made,cut,modif = IntMap.fold center_bonds ([],[],[]) (fun (made,cut,modif) i map ->
    let id1 = make_id r reid i in
    IntMap.fold map (made,cut,modif) (fun (made,cut,modif) j label ->
      let id2 = make_id r reid j in
      let s = id1 ^ "-" ^ id2 in
      match label with
        1 -> if IntSet.mem r.reactant_ids i then made, s :: cut, modif else s :: made, cut, modif
      | 2 -> made, cut, s :: modif
      | _ -> failwith "center_bonds_msg")) in
  (make_msg "BONDS MADE" made) @
  (make_msg "BONDS CUT" cut) @
  (make_msg "BONDS MODIFIED" modif)

let broken_bonds max_quality max_broken_bonds messages rl =  (* FIXME: dokończyć użycie max quality *)
(*  let time1 = Sys.time () in
  printf "broken_bonds 1 max_quality=%d max_broken_bonds=%d |rl|=%d\n%!" max_quality max_broken_bonds (Xlist.size rl);*)
(*   let max_broken_bonds = 6 in *)
  let msg = ref messages in
  let arl = Xlist.fold rl [] (fun arl r ->
    let result,history,alt_labels_list = CommonSubstructure.multilevel_label_reaction3 r [] [] in
    let alt_labels_list = if alt_labels_list = [] then [r.empty_labels] else alt_labels_list in
(*     printf "%s\n%!" (Xml.to_string_fmt (Smiles.solutions_to_xml r.record.reaction_smile [] (Xlist.map alt_labels_list (fun labels -> r, labels, [],IntMap.empty,IntMap.empty)))); *)
    (r,Collection.of_list alt_labels_list) :: arl) in
(*  let time2 = Sys.time () in
  printf "broken_bonds 2 time=%.4f |arl|=%d\n%!" (time2 -. time1) (Xlist.size arl); *)
  let no_broken_bonds, _, rgood_cands = create_good_candidates2 msg max_broken_bonds 0 arl (*r (Collection.of_list alt_labels_list)*) in
(*       let no_broken_bonds, cands, good_cands = create_good_candidates msg max_broken_bonds 0 r (Collection.of_list alt_labels_list) in *)
(*  let time3 = Sys.time () in
  printf "broken_bonds 3 time=%.4f |rgood_cands|=%d\n%!" (time3 -. time2) (Xlist.size rgood_cands);*)
  let rgood_cands = Xlist.map rgood_cands (fun (r,good_cands) ->
      let r = {r with graph=expand_hydrogens_and_fluors r.graph;
                      original_graph=expand_hydrogens_and_fluors r.original_graph;
                      reactants=Xlist.map r.reactants (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                      products=Xlist.map r.products (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                      reactant_ids=Xlist.fold r.reactants IntSet.empty (fun set m -> IntSet.union set m.ids);
                      product_ids=Xlist.fold r.products IntSet.empty (fun set m -> IntSet.union set m.ids)} in
      r,good_cands) in
(*  let time4 = Sys.time () in
  printf "broken_bonds 4 time=%.4f |rgood_cands|=%d\n%!" (time4 -. time3) (Xlist.size rgood_cands); *)
  let rminimal_good_cands, quality = select_minimal_candidates2 rgood_cands in
(*       let sorted_good_cands = sort_uniq_matchings r minimal_good_cands in *)
(*  let time5 = Sys.time () in
  printf "broken_bonds 5 time=%.4f |rminimal_good_cands|=%d\n%!" (time5 -. time4) (Xlist.size rminimal_good_cands); *)
  let rdisamb_good_cands = Xlist.map rminimal_good_cands (fun (r,cands) -> r,CommonSubstructure.disambiguate_matchings r cands) in
      let messages = !msg @ [
        (*Printf.sprintf "no broken bonds = %d" no_broken_bonds;
        Printf.sprintf "no candidates = %d" (Collection.size cands);
        Printf.sprintf "no good candidates = %d" (Collection.size good_cands);
        Printf.sprintf "no minimal good candidates = %d" (Xlist.size minimal_good_cands);
        Printf.sprintf "no disamb good candidates = %d" (Xlist.size disamb_good_cands)*)] @
        if Xlist.size rdisamb_good_cands = 0 then [(*"No matching found." Atom labelling presented.*)] else [] in
(*  let time6 = Sys.time () in
  printf "broken_bonds 6 time=%.4f |rdisamb_good_cands|=%d\n%!" (time6 -. time5) (Xlist.size rdisamb_good_cands); *)
      let solutions =
        if Xlist.size rdisamb_good_cands > 0 then
          List.flatten (Xlist.map rdisamb_good_cands (fun (r,cands) ->
          Xlist.map cands (fun (cand,matching,hydrogens) ->
(*          Printf.printf "bb \n  cand=%s \n  matching=%s \n  hydrogens=%s\n%!" (string_of_col cand) (string_of_col matching)
           (String.concat ";" (Xlist.map (IntSet.to_list hydrogens) string_of_int));*)
          let center_bonds,center_atoms = find_reaction_center_bonds2 r (Collection.to_list matching) in
          let center_bonds_modif,_ = find_reaction_center_bonds r (Collection.to_list matching) in
          let l = Xlist.fold r.reactants [] (fun l m -> set_invisible_hydrogens center_atoms l m.tree) in
          let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens center_atoms l m.tree) in
          let labels = Labels.add_invisible_list r.empty_labels l in
          let reid = make_reid r.graph (Collection.to_list matching) in
          let reid = IntSet.fold hydrogens reid IntMap.remove in
          let reid = simplify_reid labels reid in
          let reid = IntSet.fold hydrogens reid (fun reid v -> IntMap.add reid v 0) in
          let r = {r with msg = (center_bonds_msg r center_bonds reid) @ r.msg} in
          r,labels,center_bonds_modif,center_bonds,reid)))
        else [](*Xlist.map alt_labels_list (fun labels ->
          let l = Xlist.fold r.reactants [] (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let labels = Labels.add_invisible_list labels l in
          r,labels,IntMap.empty,IntMap.empty,IntMap.empty)*) in(* FIXME*)
(*  let time7 = Sys.time () in
  printf "broken_bonds 7 time=%.4f |solutions|=%d\n%!" (time7 -. time6) (Xlist.size solutions);  *)
      let solutions = ReactionClasses.disambiguate solutions in
      let solutions = ReactionClasses.select solutions in
      let solutions = ReactionClasses.remove_cycles solutions in
      let solutions = Xlist.rev_map solutions remove_stoi in
(*   let time8 = Sys.time () in
   printf "broken_bonds 8 time=%.4f |solutions|=%d\n%!" (time8 -. time7) (Xlist.size solutions);  *)
      messages, solutions, quality, no_broken_bonds

let cand_of_matching graph matching =
  let selection = Xlist.fold matching IntSet.empty (fun selection (v,w) ->
    IntSet.add (IntSet.add selection v) w) in
  Collection.of_list (Int.fold 1 (Array.length graph - 1) [] (fun cand i ->
    if IntSet.mem selection i then
      Xlist.fold (snd graph.(i)) cand (fun cand (_,q) ->
        if IntSet.mem selection q.id then (i,q.id) :: cand else cand)
    else cand))

let print_molecule m =
  Printf.printf "%s %s {%s} {%s}\n" m.smiles
    (Smiles.string_of_tree_std_id (Smiles.smile_tree_of_atom_tree m.tree))
    (String.concat "," (Xlist.rev_map (IntSet.to_list m.ids) string_of_int))
    (String.concat "," (Xlist.rev_map (IntSet.to_list m.hf_ids) string_of_int))

let print_reaction r =
  Printf.printf "%s\nReactants:\n" r.record.reaction_smile;
  Xlist.iter r.reactants print_molecule;
  (* print_graph p.reactant_graph; *)
  Printf.printf "Products:\n";
  Xlist.iter r.products print_molecule;
  (* print_graph p.product_graph; *)
  Xlist.iter r.msg print_endline;
  ()

let print_matching l =
  let l = Xlist.sort (Collection.to_list l) (fun x y -> compare (fst x) (fst y)) in
  Printf.printf "r";
  Xlist.iter l (fun (v,_) -> Printf.printf " %3d" v);
  Printf.printf "\np";
  Xlist.iter l (fun (_,v) -> Printf.printf " %3d" v);
  Printf.printf "\n"


let broken_bonds_pat messages rl =
  (* print_endline "broken_bonds_pat"; *)
  (* let msg = ref messages in *)
  (* try *)
  let rminimal_good_cands = Collection.map rl (fun (r,matching) ->
    (* print_endline "broken_bonds_pat 1a";
    print_matching (Collection.of_list matching); *)
    let ll = Xlist.map matching (fun (v,w) -> [v;w]) in
    let labels = Labels.add_list_of_lists r.empty_labels ll in
    let cand = cand_of_matching r.graph matching in
    let r2 =  {r with graph = remove_cand_pat r.graph cand} in
    let matchings = CommonSubstructure.match_reaction_full_key 0 r2 labels in
    (* print_endline "broken_bonds_pat 1b";
    print_reaction r;
    Xlist.iter (Collection.to_list matchings) print_matching;
    print_endline "broken_bonds_pat 1c"; *)
    let cands:((int * int) Collection.collection * (int * int) Collection.collection * IntSet.t) Collection.collection =
      Collection.flatten_map matchings (fun matching ->
        match_hydrogens_and_fluors cand r2 (Collection.to_list matching)) in
    (* print_endline "broken_bonds_pat 1d"; *)
    let r = {r with graph=expand_hydrogens_and_fluors r.graph;
                    original_graph=expand_hydrogens_and_fluors r.original_graph;
                    reactants=Xlist.map r.reactants (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                    products=Xlist.map r.products (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                    reactant_ids=Xlist.fold r.reactants IntSet.empty (fun set m -> IntSet.union set m.ids);
                    product_ids=Xlist.fold r.products IntSet.empty (fun set m -> IntSet.union set m.ids)} in
    (* let rminimal_good_cands, quality = select_minimal_candidates2 [r,cands] in *)
    (*Collection.of_list rminimal_good_cands*)r,cands) in
  (* with Timeout -> raise Timeout
  | e -> (msg := !msg @ [Printexc.to_string e]; Collection.empty))) in *)
  (* print_endline "broken_bonds_pat 2"; *)
  let rdisamb_good_cands = Xlist.map (Collection.to_list rminimal_good_cands) (fun (r,cands) -> r,CommonSubstructure.disambiguate_matchings r (Collection.to_list cands)) in
  (* print_endline "broken_bonds_pat 3"; *)
  let solutions =
    if Xlist.size rdisamb_good_cands > 0 then
      List.flatten (Xlist.map rdisamb_good_cands (fun (r,cands) ->
        (* print_endline "broken_bonds_pat 4";
        print_reaction r; *)
        Xlist.map cands (fun (cand,matching,hydrogens) ->
          (* print_matching matching;
          print_endline (String.concat ";" (Xlist.map (IntSet.to_list hydrogens) string_of_int)); *)
(*          Printf.printf "bb \n  cand=%s \n  matching=%s \n  hydrogens=%s\n%!" (string_of_col cand) (string_of_col matching)
           (String.concat ";" (Xlist.map (IntSet.to_list hydrogens) string_of_int));*)
          let center_bonds,center_atoms = find_reaction_center_bonds2 r (Collection.to_list matching) in
          let center_bonds_modif,_ = find_reaction_center_bonds r (Collection.to_list matching) in
          let l = Xlist.fold r.reactants [] (fun l m -> set_invisible_hydrogens center_atoms l m.tree) in
          let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens center_atoms l m.tree) in
          (* print_endline (String.concat ";" (Xlist.map (Xlist.sort l compare) string_of_int)); *)
          let labels = Labels.add_invisible_list r.empty_labels l in
          let reid = make_reid r.graph (Collection.to_list matching) in
          let reid = IntSet.fold hydrogens reid IntMap.remove in
          let reid = simplify_reid labels reid in
          let reid = IntSet.fold hydrogens reid (fun reid v -> IntMap.add reid v 0) in
          let r = {r with msg = (center_bonds_msg r center_bonds reid) @ r.msg} in
          r,labels,center_bonds_modif,center_bonds,reid)))
    else [](*Xlist.map alt_labels_list (fun labels ->
          let l = Xlist.fold r.reactants [] (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let labels = Labels.add_invisible_list labels l in
          r,labels,IntMap.empty,IntMap.empty,IntMap.empty)*) in(* FIXME*)
(*  let time7 = Sys.time () in
  printf "broken_bonds 7 time=%.4f |solutions|=%d\n%!" (time7 -. time6) (Xlist.size solutions);  *)
      let solutions = ReactionClasses.remove_cycles solutions in
      let solutions = ReactionClasses.disambiguate solutions in
      (*let solutions = ReactionClasses.select solutions in*)
      (*let solutions = Xlist.rev_map solutions remove_stoi in *)
(*   let time8 = Sys.time () in
   printf "broken_bonds 8 time=%.4f |solutions|=%d\n%!" (time8 -. time7) (Xlist.size solutions);  *)
      messages, solutions(* , quality, no_broken_bonds *)


let calculate_broken_bonds map =
  IntMap.fold map 0 (fun i _ map2 ->
    IntMap.fold map2 i (fun i _ -> function
        1 -> i+2
      | 2 -> i+1
      | _ -> failwith "calculate_broken_bonds"))


let select_minimal solutions =
  snd (Xlist.fold solutions (max_int,[]) (fun (n,found) (r,labels,center_bonds_modif,center_bonds,reid) ->
    let new_n = calculate_broken_bonds center_bonds_modif + r.broken_bonds in
    if new_n > n then n,found else
    if new_n < n then new_n,[r,labels,center_bonds_modif,center_bonds,reid] else
    n,(r,labels,center_bonds_modif,center_bonds,reid) :: found))
