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


type 'a collection = {c: 'a list; sort: bool; uniq: bool}

let of_list collection = {c=collection; sort=false; uniq=false}
let to_list collection = collection.c

let empty = {c=[]; sort=true; uniq=true}
let size collection = Xlist.size collection.c

let sum a b = {c=a.c @ b.c; sort=false; uniq=false}

let rec partition_rec selected other selector = function
    [] -> List.rev selected, List.rev other
  | elem :: collection ->
      if selector elem then partition_rec (elem :: selected) other selector collection
      else partition_rec selected (elem :: other) selector collection

let partition collection selector =
  let selected,other = partition_rec [] [] selector collection.c in
  {c=selected; sort=collection.sort; uniq=collection.uniq},
  {c=other; sort=collection.sort; uniq=collection.uniq}

let rec quotient_rec relation found = function
    {c=[]} -> List.rev found
  | {c=elem :: collection; sort=sort; uniq=uniq} ->
      let {c=equiv; sort=sort; uniq=uniq},other = partition {c=collection; sort=sort; uniq=uniq} (relation elem) in
      quotient_rec relation ({c=elem :: equiv; sort=sort; uniq=uniq} :: found) other

let quotient collection relation =
  {c=quotient_rec relation [] collection; sort=false; uniq=true}

let elem collection =
  if collection.c = [] then failwith "Collection.elem" else List.hd collection.c

let singleton elem = {c=[elem]; sort=true; uniq=true}

let is_empty collection = collection.c = []

let reverse_projection collection =
  {c=List.rev (Xlist.rev_map collection.c (fun x -> List.hd x.c)); sort=false; uniq=false}

let annotate collection annotator =
  {c=List.rev (Xlist.rev_map collection.c (fun elem ->
    elem, annotator elem)); sort=false; uniq=collection.uniq}

let deannotate collection =
  {c=List.rev (Xlist.rev_map collection.c fst); sort=false; uniq=false}

let rec combine_rec relation comp rev = function
    [],[] -> comp
  | [],_ -> failwith "combine_rec"
  | _,[] -> raise Not_found
  | elem1 :: collection1, elem2 :: collection2 ->
      if relation elem1 elem2 then combine_rec relation ((elem1,elem2) :: comp) [] (collection1, List.rev rev @ collection2) else
      combine_rec relation comp (elem2 :: rev) (elem1 :: collection1, collection2)

let combine collection1 collection2 relation =
  {c=combine_rec relation [] [] (collection1.c,collection2.c); sort=false; uniq=false}

open Xstd

let string_pair_bucket_group collection1 collection2 key_fun =
  let map = Xlist.fold collection1.c StringMap.empty (fun map x ->
    try StringMap.add_inc map (key_fun x) ([x],[]) (fun (e1,e2) -> x :: e1, e2)
    with Not_found -> map) in
  let map = Xlist.fold collection2.c map (fun map x ->
    try StringMap.add_inc map (key_fun x) ([],[x]) (fun (e1,e2) -> e1, x :: e2)
    with Not_found -> map) in
  {c=StringMap.fold map [] (fun l _ (v1,v2) ->
    ({c=List.rev v1; sort=collection1.sort; uniq=collection1.uniq},
     {c=List.rev v2; sort=collection2.sort; uniq=collection2.uniq}) :: l); sort=false; uniq=true}

let string_pair_bucket_group2 collection1 collection2 key_fun1 key_fun2 =
  let map = Xlist.fold collection1.c StringMap.empty (fun map x ->
    try StringMap.add_inc map (key_fun1 x) ([x],[]) (fun (e1,e2) -> x :: e1, e2)
    with Not_found -> map) in
  let map = Xlist.fold collection2.c map (fun map x ->
    try StringMap.add_inc map (key_fun2 x) ([],[x]) (fun (e1,e2) -> e1, x :: e2)
    with Not_found -> map) in
  {c=StringMap.fold map [] (fun l _ (v1,v2) ->
    ({c=List.rev v1; sort=collection1.sort; uniq=collection1.uniq},
     {c=List.rev v2; sort=collection2.sort; uniq=collection2.uniq}) :: l); sort=false; uniq=true}

let int_pair_bucket_group collection1 collection2 key_fun =
  let map = Xlist.fold collection1.c IntMap.empty (fun map x ->
    try IntMap.add_inc map (key_fun x) ([x],[]) (fun (e1,e2) -> x :: e1, e2)
    with Not_found -> map) in
  let map = Xlist.fold collection2.c map (fun map x ->
    try IntMap.add_inc map (key_fun x) ([],[x]) (fun (e1,e2) -> e1, x :: e2)
    with Not_found -> map) in
  {c=IntMap.fold map [] (fun l _ (v1,v2) ->
    ({c=List.rev v1; sort=collection1.sort; uniq=collection1.uniq},
     {c=List.rev v2; sort=collection2.sort; uniq=collection2.uniq}) :: l); sort=false; uniq=true}

let rec for_all collection condition =
  if collection.c = [] then true else
  if condition (List.hd collection.c) then for_all {collection with c=List.tl collection.c} condition else false

let rec exists collection condition =
  if collection.c = [] then false else
  if condition (List.hd collection.c) then true else exists {collection with c=List.tl collection.c} condition

let flatten_map collection f =
(*   print_endline "Collection.flatten_map"; *)
(*   {c=List.flatten (Xlist.rev_map collection.c (fun x -> (f x).c)); sort=false; uniq=false} *)
  {c=Xlist.fold collection.c [] (fun c x -> (f x).c @ c); sort=false; uniq=false}

let map collection f = {c=Xlist.rev_map collection.c f; sort=false; uniq=false}

let flatten collection =
  {c=List.flatten (Xlist.rev_map collection.c (fun x -> x.c)); sort=false; uniq=false}

open Big_int_Z

let int_log n = truncate (ceil (log10 (float_of_big_int n)))

let newton_symbol n k =
  let k = if k > n/2 then n-k else k in
  let x = Int.fold 1 k unit_big_int (fun x i ->
    mult_int_big_int (n-i+1) x) in
  let y = Int.fold 1 k unit_big_int (fun y i ->
    mult_int_big_int i y) in
  div_big_int x y

exception Empty

let number_of_matching_candidates ll =
  let b = Xlist.fold ll.c true (fun b (w,v) ->
    if size w * size v > 0 then false else b) in
  if b then raise Empty else
  Xlist.fold ll.c unit_big_int (fun n (w,v) ->
    let x = min (size w) (size v) in
    let y = max (size w) (size v) in
    Int.fold 0 (x-1) n (fun n k ->
      mult_int_big_int (y-k) n))

let number_of_exclusion_schemes ll =
  let b = Xlist.fold ll.c true (fun b (w,v) ->
    if size w * size v > 0 then false else b) in
  if b then raise Empty else
  Xlist.fold ll.c unit_big_int (fun n (w,v) ->
    let x = min (size w) (size v) in
    let y = max (size w) (size v) in
    mult_big_int (newton_symbol y x) n)

let rec generate_combinations_rec n l =
  if n = 0 then [[]] else
  if n = Xlist.size l then [l] else
  match l with
    [] -> failwith "generate_combinations_rec"
  | x :: l ->
    let l1 = generate_combinations_rec n l in
    let l2 = generate_combinations_rec (n-1) l in
    l1 @ (Xlist.map l2 (fun l -> x :: l))

let generate_combinations n collection =
  {c=Xlist.rev_map (generate_combinations_rec n collection.c) (fun l ->
    {c=List.rev l; sort=collection.sort; uniq=collection.uniq}); sort=false; uniq=collection.uniq}

let rec generate_partition_combinations_rec n l =
  if n = 0 then [[],l] else
  if n = Xlist.size l then [l,[]] else
  match l with
    [] -> failwith "generate_partition_combinations_rec"
  | x :: l ->
    let l1 = generate_partition_combinations_rec n l in
    let l2 = generate_partition_combinations_rec (n-1) l in
    (Xlist.rev_map l1 (fun (a,b) -> a, x :: b)) @
    (Xlist.rev_map l2 (fun (a,b) -> x :: a, b))

let generate_partition_combinations n collection =
  {c=Xlist.rev_map (generate_partition_combinations_rec n collection.c) (fun (l1,l2) ->
    {c=List.rev l1; sort=collection.sort; uniq=collection.uniq},
    {c=List.rev l2; sort=collection.sort; uniq=collection.uniq}); sort=false; uniq=collection.uniq}

let rec generate_permutations_rec = function
    [],[] -> [[]]
  | a :: la, lb -> generate_permutations_list [] a la lb
  | _ -> failwith "generate_permutations_rec"

and generate_permutations_list rev a la = function
    b :: lb -> (Xlist.rev_map (generate_permutations_rec (la,rev @ lb)) (fun l -> (a,b) :: l)) @ (generate_permutations_list (b :: rev) a la lb)
  | [] -> []

let generate_permutations (collection,ets) =
  {c=Xlist.rev_map (generate_permutations_rec (collection.c,ets.c)) (fun l -> {c=l; sort=false; uniq=false}); sort=false; uniq=false}

let generate_product_permutations ll =
(*   print_endline "generate_product_permutations 1"; *)
  let l = Xlist.fold ll.c [] (fun l (la,lb) -> (generate_permutations_rec (la.c,lb.c)) :: l) in
(*   print_endline "generate_product_permutations 2"; *)
  {c=Xlist.rev_map (Xlist.multiply_list l) (fun l -> {c=List.flatten l; sort=false; uniq=false}); sort=false; uniq=false}

let rec factorial x n =
  if n < 2 then x else factorial (mult_int_big_int n x) (n-1)

let number_of_product_permutations ll =
  Xlist.fold ll.c unit_big_int (fun n (la,lb) ->
    mult_big_int n (factorial unit_big_int (size la)))

(*let rec select_one found rev = function
    [] -> found
  | x :: l -> select_one ((x, rev @ l) :: found) (x :: rev) l

let rec generate_permutations =
    [] -> [[]]
  | l -> Xlist.fold (select_one [] [] l) [] (fun perms (first,rest) ->
           Xlist.fold (generate_permutations rest) perms (fun perms perm ->
             (first :: perm) :: perms))

(* l - lista kolekcji obiektów
   ll - lista permutacji *)
let generate_product_permutations l =
  let ll = Xlist.map l generate_permutations in
  Xlist.map (Xlist.multiply_list ll) List.flatten   *)

let rec number_of_products_rec = function
    [] -> unit_big_int
  | l :: ll ->
     let x = number_of_products_rec ll in
     mult_int_big_int (size l) x

let number_of_products col =
  number_of_products_rec col.c

let product col =
  let ll = (*List.rev*) (Xlist.rev_map col.c to_list) in
  {c=Xlist.rev_map (Xlist.multiply_list ll) (fun l -> {c=l; sort=false; uniq=false}); sort=false; uniq=false}

let flatten_product col =
  let ll = (*List.rev*) (Xlist.rev_map col.c to_list) in
  {c=Xlist.rev_map (Xlist.multiply_list ll) (fun l -> {c=List.flatten (Xlist.rev_map l to_list); sort=false; uniq=false}); sort=false; uniq=false}

let get_exclusion_schemas_complete ll =
  let l1,l2 = Xlist.fold ll.c ([],[]) (fun (l1,l2) (w,v) ->
    if size w = size v then l1,l2 else
    if size w < size v then
     l1, (generate_combinations_rec (size v - size w) v.c) :: l2
    else (generate_combinations_rec (size w - size v) w.c) :: l1, l2) in
  let l1 = Xlist.map (Xlist.multiply_list l1) (fun l -> {c=List.flatten l; sort=false; uniq=false}) in
  let l2 = Xlist.map (Xlist.multiply_list l2) (fun l -> {c=List.flatten l; sort=false; uniq=false}) in
(*   (if l1 = [] then [[]] else l1), (if l2 = [] then [[]] else l2) *)
  {c=l1; sort=false; uniq=false}, {c=l2; sort=false; uniq=false}

let get_exclusion_schemas ll =
  let l = Xlist.fold ll.c [] (fun l (w,v) ->
    if size w = size v then l else
    if size w < size v then
     (generate_combinations_rec (size v - size w) v.c) :: l
    else (generate_combinations_rec (size w - size v) w.c) :: l) in
  {c=Xlist.map (Xlist.multiply_list l) (fun l -> {c=List.flatten l; sort=false; uniq=false}); sort=false; uniq=false}

let rec get_border_schemas_rec matching found = function
    [] -> {c=matching; sort=false; uniq=false} :: found
  | (w,v) :: l ->
    let lw,lv = if size w < size v then
     [w.c], generate_combinations_rec (size w) v.c
    else generate_combinations_rec (size v) w.c, [v.c] in
    Xlist.fold lw found (fun found w ->
      Xlist.fold lv found (fun found v ->
        get_border_schemas_rec (({c=w; sort=false; uniq=false},{c=v; sort=false; uniq=false}) :: matching) found l))

let get_border_schemas ll matching =
  {c=get_border_schemas_rec matching.c [] ll.c; sort=false; uniq=false}

let sort collection = {c=Xlist.sort collection.c compare; sort=true; uniq=collection.uniq}

let rec uniq_rec = function
    x :: y :: l -> if x = y then uniq_rec (x :: l) else x :: (uniq_rec (y :: l))
  | l -> l

let uniq collection =
  let collection = if collection.sort then collection else sort collection in
  {c=uniq_rec collection.c; sort=true; uniq=true}

let rec inter rev = function
    [],_ -> List.rev rev
  | _,[] -> List.rev rev
  | x1 :: l1, x2 :: l2 ->
      if x1 = x2 then inter (x1 :: rev) (l1,l2) else
      if x1 < x2 then inter rev (l1,x2 :: l2) else
      inter rev (x1 :: l1,l2)

let intersection collection1 collection2 =
(*   print_endline "Collection.intersection1"; *)
  let collection1 = if collection1.sort then collection1 else sort collection1 in
(*   print_endline "Collection.intersection2"; *)
  let collection2 = if collection2.sort then collection2 else sort collection2 in
(*   print_endline "Collection.intersection3"; *)
  {c=inter [] (collection1.c,collection2.c); sort=true; uniq=collection1.uniq || collection2.uniq}

let rec group_rec rev = function
    (k1,v1) :: (k2,v2) :: l -> if k1 = k2 then group_rec rev ((k1,v2@v1) :: l) else group_rec ((k1,v1) :: rev) ((k2,v2) :: l)
  | l -> List.rev (List.rev l @ rev)

let group collection =
  let collection = sort {c=Xlist.rev_map collection.c (fun (k,v) -> k,[v]); sort=false; uniq=collection.uniq} in
  {c=Xlist.rev_map (group_rec [] collection.c) (fun (k,v) -> k,{c=v; sort=false; uniq=false}); sort=false; uniq=true}

let group2 collection =
  let collection = sort {c=Xlist.rev_map collection.c (fun (v,k) -> k,[v]); sort=false; uniq=collection.uniq} in
  {c=Xlist.rev_map (group_rec [] collection.c) (fun (k,v) -> {c=v; sort=false; uniq=false},k); sort=false; uniq=true}

let to_string collection delim elt_fun =
  String.concat delim (Xlist.map collection.c elt_fun)

let to_string_as_list collection elt_fun =
  "[" ^ to_string collection ";" elt_fun ^ "]"

let print_cc collection name name2 elt_fun =
  print_endline name;
  if name2 = "" then
    Xlist.iter collection.c (fun col ->
      print_endline (to_string col " " elt_fun))
  else
    Xlist.iter collection.c (fun col ->
      print_endline ("  " ^ name2);
      print_endline ("  " ^ to_string col "\n  " elt_fun))

let print_ccc collection name name2 name3 elt_fun =
  print_endline name;
  Xlist.iter collection.c (fun col ->
    print_endline ("  " ^ name2);
    if name3 = "" then
    Xlist.iter col.c (fun c ->
      print_endline ("  " ^ to_string c " " elt_fun))
    else
      Xlist.iter col.c (fun c ->
        print_endline ("    " ^ name3);
        print_endline ("    " ^ to_string c "\n    " elt_fun)))
