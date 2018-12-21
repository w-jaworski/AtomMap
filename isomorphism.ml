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

(* algorytm maczowania:
   mamy zbior dopasowan S;
   dopasowanie jest zbiorem parnalezacych do jednego z 3 typow:
     - para wierzchołków, która ma dopasowanych wszystkich sasiadów
     - para wierzchołków, która ma dopasowaną część sasiadów i nie została wykonana próba dopasowania pozostałych
     - para wierzchołków, która mają niedopasowywalnych sąsiadów
   wybierz dopasowanie D z S, które ma parę wierzchołków P z niedopasowywanymi jeszcze sąsiadami
   dopasuj jednego z sąsiadów P1 na wszystkie możliwe sposoby
     - dodaj do S dopasowanie D rozszerzone o tego sąsiada na wszystkie sposoby
     - jeśli nie ma żadnego sposobu uznaj go za niedopasowywalnego

   dopasowanie polega na
     - sprawdzeniu etykiet wierzchołków
     - sprawdzeniu czy wierzchołki nie są już dopasowane
*)

let translate_graph bond_flag hash graph selection =
  let selection = IntSet.of_list (Collection.to_list selection) in
  (* Printf.printf "|selection|=%d |graph|=%d\n" (IntSet.size selection) (Array.length graph); *)
  let map = Array.fold_left (fun map (p,l) ->
    if IntSet.mem selection p.id then
      IntMap.add map p.id (hash.(p.id),l)
    else map) IntMap.empty graph in
  IntMap.map map (fun (s,l) ->
    s, Xlist.fold l [] (fun l (b,p) ->
      if IntSet.mem selection p.id then
        ((if bond_flag then Hash.key_of_bond b else "") ^ fst (IntMap.find map p.id), p.id) :: l
      else l))

(****
let select_matching_root_single hash selection =
  let map = IntSet.fold selection StringQMap.empty (fun map n ->
    StringQMap.add map hash.(n)) in
  fst (StringQMap.fold map ("",1000000) (fun (x,v) s n ->
    if n < v then s,n else x,v))
****)
let select_matching_root hash1 hash2 selection1 selection2 =
(*   let ll = Collection.to_list (Collection.string_pair_bucket_group selection1 selection2 (fun n -> hash1.(n)) (fun n -> hash2.(n))) *)
  let map = Hash.hash_bucket_group2 hash1 hash2 selection1 selection2 in
  fst (StringMap.fold map ("",1000000) (fun (x,v) s (l1,l2) ->
    let v1 = Xlist.size l1 in
    let v2 = Xlist.size l2 in
    if v1 * v2 = 0 then x,v else
    if v1 * v2 < v then s, v1 * v2 else x,v))
(****
let select_multiple_matching_roots conn hash1 hash2 selection1 selection2 =
  let map = Hash.hash_bucket_group2 hash1 hash2 selection1 selection2 in
  Xlist.map conn (fun l ->
    let set = IntSet.of_list l in
    fst (StringMap.fold map ("",1000000) (fun (x,v) s (l1,l2) ->
      let set1 = IntSet.of_list l1 in
      if IntSet.size (IntSet.intersection set set1) = 0 then x,v else
      let v1 = Xlist.size l1 in
      let v2 = Xlist.size l2 in
      if v1 * v2 = 0 then x,v else
      if v1 * v2 < v then s, v1 * v2 else x,v)),
    Xlist.size l)

****)
let create_starting_assumptions_exact s ma mb =
   (* print_endline s; *)
  let la = IntMap.fold ma [] (fun la n (r,_) -> if s = r then n :: la else la) in
  let lb = IntMap.fold mb [] (fun lb n (r,_) -> if s = r then n :: lb else lb) in
  (* print_endline (string_of_selection la); *)
  (* print_endline (string_of_selection lb); *)
  if Xlist.size la < Xlist.size lb then (* zakładam, że wszyskie wyróżnione atomy z cząsteczki, w której jest ich mniej przechodzą do drugiej cząsteczki *)
    Xlist.map lb (fun b -> [List.hd la, b])
  else Xlist.map la (fun a -> [a, List.hd lb])
(****
let create_starting_assumptions_exact2 s ma mb =
  let la = IntMap.fold ma [] (fun la n (r,_) -> if s = r then n :: la else la) in
  let lb = IntMap.fold mb [] (fun lb n (r,_) -> if s = r then n :: lb else lb) in
  match la,lb with
    [a1],[b1] -> [[a1,b1]]
  | [a1;a2],[b1;b2] -> [[a1,b1;a2,b2];[a1,b2;a2,b1]]
  | _ -> failwith "create_starting_assumptions_exact2"

(*let rec random_select_rec selected x =
  let i = Random.int x in
  if IntSet.mem selected i then random_select_rec selected x else i

let random_select n l =
  let a = Array.of_list l in
  let selected = Int.fold 1 n IntSet.empty (fun selected _ ->
    IntSet.add selected (random_select_rec selected (Xlist.size l))) in
  IntSet.fold selected [] (fun l i -> a.(i) :: l)

let rec make_random_pairs_exact n = function
    [],[] -> [[]]
  | a :: la, lb -> random_select n (make_random_pairs_exact_list n [] a la lb)
  | _ -> failwith "make_random_pairs_exact"

and make_random_pairs_exact_list n rev a la = function
    b :: lb -> (Xlist.map (make_random_pairs_exact n (la,rev @ lb)) (fun l -> (a,b) :: l)) @ (make_random_pairs_exact_list n (b :: rev) a la lb)
  | [] -> []*)

let rec random_permutation selected x n =
(*   Printf.printf "%s x=%d n=%d\n" (string_of_int_list (IntSet.to_list selected)) x n; *)
  if n = 0 then [] else
  let i = Random.int x in
  if IntSet.mem selected i then random_permutation selected x n else
  i :: (random_permutation (IntSet.add selected i) x (n-1))

(*let rec random_permutation_rec n =
  if n = 0 then [] else (Random.int n) :: (random_permutation_rec (n-1))

let rec random_permutation2 = function
    [] -> []
  | n :: l -> n :: (random_permutation2 (Xlist.map l (fun i -> if i >= n then i+1 else i)))

let random_permutation n =
  let l = random_permutation_rec n in
  Printf.printf "a %s\n" (string_of_int_list l);
  random_permutation2 l*)

let make_random_pairs_exact n (l,l2) =
(*   Printf.printf "make_random_pairs_exact %d %s %s\n%!" n (string_of_int_list l) (string_of_int_list l2) ; *)
  let a = Array.of_list l in
(*   let a2 = Array.of_list l2 in *)
  Int.fold 1 n [] (fun found _ ->
    let p = random_permutation IntSet.empty (Xlist.size l) (Xlist.size l) in
(*     Printf.printf "b %s\n" (string_of_int_list p);  *)
    (Xlist.fold2 p l2 [] (fun sol i w -> (a.(i),w) :: sol)) :: found)
****)
let get_possible_matchings_exact ma mb va vb matcheda matchedb =
  if fst (IntMap.find ma va) = "R" then Collection.singleton Collection.empty else (
  (* print_endline "get_possible_matchings_exact 1"; *)
  let la = snd (IntMap.find ma va) in
  let lb = snd (IntMap.find mb vb) in
  (* Printf.printf "la=%s\n" (String.concat " " (Xlist.map la (fun (s,n) -> s ^ "-" ^ string_of_int n))); *)
  let la,lra = Xlist.fold la ([],[]) (fun (la,lr) (s,n) -> if s = "-R" then la, n :: lr else (s,n) :: la, lr) in
  (* Printf.printf "lra=%s\n" (string_of_selection lra); *)
  let seta = Xlist.fold la IntSet.empty (fun set (s,n) -> if IntMap.mem matcheda n then IntSet.add set n else set) in
  let setb = Xlist.fold lb IntSet.empty (fun set (s,n) -> if IntMap.mem matchedb n then IntSet.add set n else set) in
  let ba = IntSet.fold seta true (fun b n -> (IntSet.mem setb (IntMap.find matcheda n)) && b) in
  let bb = IntSet.fold setb true (fun b n -> (IntSet.mem seta (IntMap.find matchedb n)) && b) in
  if (not ba) || (not bb) then Collection.empty else (
  (* print_endline "get_possible_matchings_exact 2"; *)
  let la = Xlist.fold la [] (fun l (s,n) -> if IntSet.mem seta n then l else (n,s) :: l) in
  let lb = Xlist.fold lb [] (fun l (s,n) -> if IntSet.mem setb n then l else (n,s) :: l) in
  if Xlist.size lra > Xlist.size lb then Collection.empty else (
  (* print_endline "get_possible_matchings_exact 2"; *)
  (* Printf.printf "lra=%s\n" (string_of_selection lra); *)
  let cands = Collection.generate_combinations (Xlist.size lra) (Collection.of_list lb) in
  Collection.flatten_map cands (fun lrb ->
    let lrb = Xlist.rev_map (Collection.to_list lrb) fst in
    (* Printf.printf "lrb=%s\n" (string_of_selection lrb); *)
    let lb = Xlist.fold lb [] (fun l (n,s) -> if Xlist.mem lrb n then l else (n,s) :: l) in
    let ll = Collection.string_pair_bucket_group (Collection.of_list la) (Collection.of_list lb) snd in
    let ll = Collection.of_list ((Collection.of_list lra,Collection.of_list lrb) :: Xlist.map (Collection.to_list ll) (fun (x,y) -> Collection.deannotate x, Collection.deannotate y)) in
    if Collection.exists ll (fun (la,lb) -> Collection.size la <> Collection.size lb) then Collection.empty else
  (*let ll =*) Collection.generate_product_permutations ll)))) (*in
  Xlist.map (Collection.to_list ll) Collection.to_list*)
(*  let map = Xlist.fold la StringMap.empty (fun map (s,n) ->
    StringMap.add_inc map s ([n],[]) (fun (l1,l2) -> n :: l1, l2)) in
  let map = Xlist.fold lb map (fun map (s,n) ->
    StringMap.add_inc map s ([],[n]) (fun (l1,l2) -> l1, n :: l2)) in
  if StringMap.fold map false (fun b _ (la,lb) -> Xlist.size la <> Xlist.size lb || b) then [] else
  let l = StringMap.fold map [] (fun l _ (la,lb) -> (make_pairs_exact (la,lb)) :: l) in
  Xlist.map (Xlist.multiply_list l) List.flatten*)
(****
let get_possible_matchings_exact_prematched (prematcheda,prematchedb) ma mb va vb matcheda matchedb =
  let la = snd (IntMap.find ma va) in
  let lb = snd (IntMap.find mb vb) in
(*   Printf.printf "g1 la=%s lb=%s\n" (string_of_int_list (Xlist.map la snd)) (string_of_int_list (Xlist.map lb snd)); *)
  let la = Xlist.fold la [] (fun l (s,n) -> if IntSet.mem prematcheda n then l else (s,n) :: l) in
  let lb = Xlist.fold lb [] (fun l (s,n) -> if IntSet.mem prematchedb n then l else (s,n) :: l) in
(*   Printf.printf "g2 la=%s lb=%s\n" (string_of_int_list (Xlist.map la snd)) (string_of_int_list (Xlist.map lb snd)); *)
  let seta = Xlist.fold la IntSet.empty (fun set (s,n) -> if IntMap.mem matcheda n then IntSet.add set n else set) in
  let setb = Xlist.fold lb IntSet.empty (fun set (s,n) -> if IntMap.mem matchedb n then IntSet.add set n else set) in
  let ba = IntSet.fold seta true (fun b n -> (IntSet.mem setb (IntMap.find matcheda n)) && b) in
  let bb = IntSet.fold setb true (fun b n -> (IntSet.mem seta (IntMap.find matchedb n)) && b) in
  if (not ba) || (not bb) then (*(Printf.printf "g1 va=%d vb=%d\n" va vb;*)[] else
  let la = Xlist.fold la [] (fun l (s,n) -> if IntSet.mem seta n then l else (s,n) :: l) in
  let lb = Xlist.fold lb [] (fun l (s,n) -> if IntSet.mem setb n then l else (s,n) :: l) in
(*   Printf.printf "g3 la=%s lb=%s\n" (string_of_int_list (Xlist.map la snd)) (string_of_int_list (Xlist.map la snd)); *)
  let ll = Pair.string_bucket_group (Pair.singleton la lb) fst snd in
  if Pair.fold ll false (fun b (la,lb) -> Xlist.size la <> Xlist.size lb || b) then [] else
  Pair.make_matchings ll
(*  let map = Xlist.fold la StringMap.empty (fun map (s,n) ->
    StringMap.add_inc map s ([n],[]) (fun (l1,l2) -> n :: l1, l2)) in
  let map = Xlist.fold lb map (fun map (s,n) ->
    StringMap.add_inc map s ([],[n]) (fun (l1,l2) -> l1, n :: l2)) in
  if StringMap.fold map false (fun b _ (la,lb) -> Xlist.size la <> Xlist.size lb || b) then (*(Printf.printf "g2 va=%d vb=%d\n" va vb;*)[] else
  let l = StringMap.fold map [] (fun l _ (la,lb) -> (make_pairs_exact (la,lb)) :: l) in
  Xlist.map (Xlist.multiply_list l) List.flatten*)
****)
type decision_tree =
    T of int * int * decision_tree Collection.collection
  | Tfinish
  | Terror

let counter = ref 0
let max_tree_size = 1000000
exception TreeToBig

let rec create_decision_tree matching_fun matching_size ma mb matcheda matchedb = function
    [] -> incr counter; if matching_size = 0 then Tfinish else Terror
  | (va,vb) :: vertices ->
       (* Printf.printf "va=%d vb=%d %d matching_size=%d\n" va vb (Xlist.size vertices) matching_size; *)
       if Sys.time () -. !Types.time > Types.timeout then ((*print_endline "TIMEOUT";*) raise Timeout) else
      incr counter;
      if !counter > max_tree_size then raise TreeToBig;
      if matching_size = 0 then Tfinish else
      if IntMap.mem matcheda va then
        if IntMap.find matcheda va = vb then create_decision_tree matching_fun matching_size ma mb matcheda matchedb vertices else Terror else
      if IntMap.mem matchedb vb then
        if IntMap.find matchedb vb = va then create_decision_tree matching_fun matching_size ma mb matcheda matchedb vertices else Terror else
      let matcheda = IntMap.add matcheda va vb in
      let matchedb = IntMap.add matchedb vb va in
      let ll = matching_fun ma mb va vb matcheda matchedb in
      if Collection.is_empty ll && matching_size > 1 then Terror else
      T(va,vb,Collection.map ll (fun l ->
        create_decision_tree matching_fun (matching_size - 1) ma mb matcheda matchedb ((Collection.to_list l) @ vertices)))

let rec extract_matchings matched = function
    T(va,vb,l) -> Collection.flatten_map l (extract_matchings ((va,vb) :: matched))
  | Tfinish -> Collection.singleton matched
  | Terror -> Collection.empty

let rec has_matching = function
    T(_,_,l) -> Collection.exists l has_matching
  | Tfinish -> true
  | Terror -> false
(****
let find_symmetries hash g selection =
  let g = translate_graph hash.(3) g selection in
  let matching_root = select_matching_root_single hash.(3) selection in
  let l = create_starting_assumptions_exact matching_root g g in
  counter := 0;
  let trees = Xlist.map l (create_decision_tree get_possible_matchings_exact (IntMap.size g) g g IntMap.empty IntMap.empty) in
  List.flatten (Xlist.map trees (extract_matchings []))
****)
let rec are_isomorphic_rec g1 g2 = function
    [] -> (*print_endline "not isomorphic";*) false
  | x :: l ->
(*       Printf.printf "are_isomorphic_rec\n";  *)
      if has_matching (create_decision_tree get_possible_matchings_exact (IntMap.size g1) g1 g2 IntMap.empty IntMap.empty x) then
        ((*Printf.printf "are_isomorphic_rec: matching found\n";*) true) else
      are_isomorphic_rec g1 g2 l

let are_isomorphic bond_flag level hash1 hash2 g1 g2 selection1 selection2 =
(*      Printf.printf "selection1=%s\n" (string_of_selection (IntSet.to_list selection1));
      Printf.printf "selection2=%s\n" (string_of_selection (IntSet.to_list selection2));  *)
  if Collection.is_empty selection1 && Collection.is_empty selection2 then true else
  let g1 = translate_graph bond_flag hash1.(min level 1) g1 selection1 in
  let g2 = translate_graph bond_flag hash2.(min level 1) g2 selection2 in
  let matching_root = select_matching_root hash1.(min level 1) hash2.(min level 1) selection1 selection2 in
  let l = create_starting_assumptions_exact matching_root g1 g2 in
  counter := 0;
  are_isomorphic_rec g1 g2 l

let rec are_isomorphic_lists_rec bond_flag level = function
    [],[] -> true
  | _,[] -> failwith "are_isomorphic_lists"
  | [],_ -> failwith "are_isomorphic_lists"
  | (hash1,graph1,selection1) :: hgsl1, (hash2,graph2,selection2) :: hgsl2 ->
      if are_isomorphic bond_flag level hash1 hash2 graph1 graph2 selection1 selection2 then
        are_isomorphic_lists_rec bond_flag level (hgsl1,hgsl2)
      else false

let rec are_isomorphic_lists_rec2 bond_flag level graph1 graph2 hash1 hash2 = function
    [] -> true
  | (selection1, selection2) :: l ->
      if are_isomorphic bond_flag level hash1 hash2 graph1 graph2 selection1 selection2 then
        are_isomorphic_lists_rec2 bond_flag level graph1 graph2 hash1 hash2 l
      else false

let are_isomorphic_lists bond_flag level (_,hgsl1) (_,hgsl2) =
  are_isomorphic_lists_rec bond_flag level (hgsl1,hgsl2)

let rec are_isomorphic_alt_lists bond_flag level (x1,(graph1,hash1)) (x2,(graph2,hash2)) = function
    [] -> false
  | l :: alt_list -> if are_isomorphic_lists_rec2 bond_flag level graph1 graph2 hash1 hash2 (Collection.to_list l) then true else
                       are_isomorphic_alt_lists bond_flag level (x1,(graph1,hash1)) (x2,(graph2,hash2)) alt_list

let find_isomorphisms bond_flag level hash1 hash2 g1 g2 selection1 selection2 =
(*      Printf.printf "selection1=%s\n" (string_of_selection (IntSet.to_list selection1));
      Printf.printf "selection2=%s\n" (string_of_selection (IntSet.to_list selection2));  *)
  (* Printf.printf "|hash1|=%d |hash2|=%d %s %s\n%!" (Array.length hash1.(0)) (Array.length hash2.(0))
    (string_of_selection (Collection.to_list selection1))
    (string_of_selection (Collection.to_list selection2)); *)
  let g1 = translate_graph bond_flag hash1.(min level 1) g1 selection1 in
  let g2 = translate_graph bond_flag hash2.(min level 1) g2 selection2 in
  (* Printf.printf "|g1|=%d |g2|=%d\n%!" (IntMap.size g1) (IntMap.size g2); *)
      (* print_endline "eeee1"; *)
  let matching_root = select_matching_root hash1.(min level 1) hash2.(min level 1) selection1 selection2 in
      (* print_endline "eeee4"; *)
  let l = create_starting_assumptions_exact matching_root g1 g2 in
  counter := 0;
      (* print_endline "eeee2"; *)
  Collection.flatten_map (Collection.of_list l) (fun x ->
    let tree = create_decision_tree get_possible_matchings_exact (IntMap.size g1) g1 g2 IntMap.empty IntMap.empty x in
      (* print_endline "eeee3"; *)
    Collection.map (extract_matchings [] tree) Collection.of_list)
(****
let extend_matching level hash1 hash2 g1 g2 selection1 selection2 matching =
  let g1 = translate_graph hash1.(min level 1) g1 selection1 in
  let g2 = translate_graph hash2.(min level 1) g2 selection2 in
  let trees = create_decision_tree
    (get_possible_matchings_exact_prematched (Xlist.fold matching (IntSet.empty,IntSet.empty) (fun (set1,set2) (v,w) -> IntSet.add set1 v, IntSet.add set2 w)))
    (IntMap.size g1) g1 g2 IntMap.empty IntMap.empty matching in
  extract_matchings [] trees

let select_single_matching_roots hash1 hash2 selection1 selection2 =
  let map = Hash.hash_bucket_group2 hash1 hash2 selection1 selection2 in
  StringMap.fold map [] (fun l _ -> function
      [v],[w] -> (v,w) :: l
    | _ -> l)

let extend_matching2 level hash1 hash2 g1 g2 selection1 selection2 matching = failwith "ni"
(*  let e1 = Xlist.fold matching [] (fun e1 (v,_) -> v :: e1) in
  let selection1,selection2 = Xlist.fold matching (selection1,selection2) (fun (selection1,selection2) (v,w) -> IntSet.remove selection1 v, IntSet.remove selection2 w) in
(*   Printf.printf "selection1=%s\n" (string_of_int_list (IntSet.to_list selection1));  *)
  let conn = get_connected_components e1 g1 selection1 in
(*  Printf.printf "selection1=%s\n" (string_of_int_list (IntSet.to_list selection1));
  Printf.printf "%s\n" (String.concat " " (Xlist.map conn string_of_int_list));*)
  let g1 = translate_graph hash1.(min level 1) g1 selection1 in
  let g2 = translate_graph hash2.(min level 1) g2 selection2 in
  let matching_roots = select_multiple_matching_roots conn hash1.(min level 1) hash2.(min level 1) selection1 selection2 in
  let l = Xlist.map matching_roots (fun (matching_root,size) ->
    let l = create_starting_assumptions_exact(*2*) matching_root g1 g2 in  (* FIXME *)
(*     Printf.printf "size=%d %s %s\n" size matching_root (String.concat " " (Xlist.map l (fun l -> "[" ^ (String.concat " " (Xlist.map l (fun (w,v) -> Printf.sprintf "%d->%d" w v)) ^ "]"))));   *)
    let trees = Xlist.map l (fun x -> create_decision_tree get_possible_matchings_exact size g1 g2 IntMap.empty IntMap.empty x) in
    let found = List.flatten (Xlist.rev_map trees (extract_matchings [])) in
(*     Printf.printf "[%s]\n" (String.concat " " (Xlist.map found (fun matching -> string_of_int (Xlist.size matching)))); *)
    found) in
  Xlist.rev_map (Xlist.multiply_list ([matching] :: l)) List.flatten*)

let rec remove_conns set = function
    [] -> []
  | c :: conn ->
      let b = Xlist.fold c true (fun b i -> b && (not (IntSet.mem set i))) in
      if b then c :: (remove_conns set conn) else remove_conns set conn

let rec extend_matching_starting_assumptions conn2 hash matchings = function
    [] -> matchings
  | c :: conn1 ->
      let map = Xlist.fold c StringMap.empty (fun map i -> StringMap.add_inc map hash.(i) [i] (fun l -> i :: l)) in
      let x,_ = StringMap.fold map (-1,1000000) (fun (x,v) _ l -> if Xlist.size l < v then List.hd l, Xlist.size l else x, v) in
      let l = Xlist.fold conn2 [] (fun l c2 ->
        if Xlist.size c2 = Xlist.size c then
          Xlist.fold c2 l (fun l i -> if hash.(i) = hash.(x) then [x,i] :: l else l)
        else l) in
      extend_matching_starting_assumptions conn2 hash (l :: matchings) conn1

let extend_matching3 hash g selection1 selection2 conn1 conn2 division =
  let matched1,matched2 = Pair.fold division (IntSet.empty,IntSet.empty) (fun (matched1,matched2) (l1,l2) ->
    Xlist.fold l1 matched1 IntSet.add, Xlist.fold l2 matched2 IntSet.add) in
  let conn1 = remove_conns matched1 conn1 in
  let conn2 = remove_conns matched1 conn2 in
  let matchings = Pair.make_matchings division in
  let matchings = extend_matching_starting_assumptions conn2 hash [matchings] conn1 in
  let matchings = Xlist.map (Xlist.multiply_list matchings) List.flatten in
  let g1 = translate_graph hash g selection1 in
  let g2 = translate_graph hash g selection2 in
  let size = IntSet.size selection1 in
  let trees = Xlist.map matchings (fun matching ->
(*     Printf.printf "size=%d [%s]\n" size (String.concat " " (Xlist.map matching (fun (w,v) -> Printf.sprintf "%d->%d" w v))); *)
    create_decision_tree (get_possible_matchings_exact_prematched (matched1,matched2)) size g1 g2 IntMap.empty IntMap.empty matching) in
  List.flatten (Xlist.rev_map trees (extract_matchings []))
****)
