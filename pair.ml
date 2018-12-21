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


type 'a t = ('a list * 'a list) list

let empty = []
let singleton x y =  [x,y]

let to_pair_list x = x
let of_pair_list x = x

let add ll x = x :: ll
let append ll ll2 = ll @ ll2

open Xstd

let string_bucket_group ll key_fun elt_fun =
  List.flatten (Xlist.map ll (fun (l1,l2) ->
    let map = Xlist.fold l1 StringMap.empty (fun map x ->
      try StringMap.add_inc map (key_fun x) ([elt_fun x],[]) (fun (e1,e2) -> (elt_fun x) :: e1, e2)
      with Not_found -> map) in
    let map = Xlist.fold l2 map (fun map x ->
      try StringMap.add_inc map (key_fun x) ([],[elt_fun x]) (fun (e1,e2) -> e1, (elt_fun x) :: e2)
      with Not_found -> map) in
    StringMap.fold map [] (fun l _ v -> v :: l)))

let int_bucket_group ll key_fun elt_fun =
  List.flatten (Xlist.map ll (fun (l1,l2) ->
    let map = Xlist.fold l1 IntMap.empty (fun map x ->
      try IntMap.add_inc map (key_fun x) ([elt_fun x],[]) (fun (e1,e2) -> (elt_fun x) :: e1, e2)
      with Not_found -> map) in
    let map = Xlist.fold l2 map (fun map x ->
      try IntMap.add_inc map (key_fun x) ([],[elt_fun x]) (fun (e1,e2) -> e1, (elt_fun x) :: e2)
      with Not_found -> map) in
    IntMap.fold map [] (fun l _ v -> v :: l)))


let rec make_pairs_exact = function
    [],[] -> [[]]
  | a :: la, lb -> make_pairs_exact_list [] a la lb
  | _ -> failwith "make_pairs_exact"

and make_pairs_exact_list rev a la = function
    b :: lb -> (Xlist.map (make_pairs_exact (la,rev @ lb)) (fun l -> (a,b) :: l)) @ (make_pairs_exact_list (b :: rev) a la lb)
  | [] -> []

let make_matchings ll =
  let l = Xlist.fold ll [] (fun l (la,lb) -> (make_pairs_exact (la,lb)) :: l) in
  Xlist.map (Xlist.multiply_list l) List.flatten

let fold = Xlist.fold

let list_fold (l1,l2) (s1,s2) f =
  let l1 = Xlist.fold l1 s1 f in
  let l2 = Xlist.fold l2 s2 f in
  l1,l2

let list_map (l1,l2) f =
  Xlist.map l1 f, Xlist.map l2 f

(**********************************************************************************************)

open Big_int

let rec factorial n = function
    0 -> n
  | i -> factorial (mult_int_big_int i n) (i-1)


(*let string_of_int_list l =
  Printf.sprintf "[%s]" (String.concat ";" (Xlist.map l string_of_int))*)

let count_matchings ll =
  fold ll unit_big_int (fun n (*k*) (w,v) ->
    if Xlist.size w <> Xlist.size v then failwith ("count_matchings: " (*^ string_of_int_list w ^ " " ^ string_of_int_list v*)) else
    factorial n (Xlist.size w))
