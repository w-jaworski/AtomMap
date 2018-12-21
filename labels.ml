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


open Xstd

type t = int array

(* let empty = [| |] *)

let make size =
  let labels = Array.make size 0 in
  labels.(0) <- 1;
  labels

let find labels i =
  if i = 0 then failwith "Labels.find: invalid key" else
  if labels.(i) = 0 then raise Not_found else
  labels.(i)

let get labels i =
  if i = 0 then failwith "Labels.get: invalid key" else
  labels.(i)

let set labels i =
  if i = 0 then failwith "Labels.set: invalid key" else
  labels.(i) <- labels.(0)

let set_if_new labels i =
  if i = 0 then failwith "Labels.set_if_new: invalid key" else
  if labels.(i) = 0 then labels.(i) <- labels.(0) else ()

let next labels =
  labels.(0) <- labels.(0) + 1

let get_next labels =
  labels.(0)

let mem labels i =
  if i = 0 then failwith "Labels.mem: invalid key" else
  if labels.(i) = 0 then false else true

let fold labels s f =
  Int.fold 1 (Array.length labels - 1) s (fun s i ->
    if labels.(i) = 0 then s else
    f s labels.(i))

let foldi labels s f =
  Int.fold 1 (Array.length labels - 1) s (fun s i ->
    if labels.(i) = 0 then s else
    f s i labels.(i))

let add_list_of_list_pairs labels ll =
  let labels = Array.copy labels in
  Xlist.iter (Pair.to_pair_list ll) (fun (l1,l2) ->
    Xlist.iter l1 (set_if_new labels);
    Xlist.iter l2 (set_if_new labels);
    next labels);
  labels

(*let make_list_of_lists size ll =
  let labels = make size in
  Xlist.iter ll (fun l ->
    Xlist.iter l (set labels);
    next labels);
  labels*)

let add_list_of_lists labels ll =
  let labels = Array.copy labels in
  Xlist.iter ll (fun l ->
    Xlist.iter l (set labels);
    next labels);
  labels

let map labels nlabels next =
  let labels = Array.copy labels in
  labels.(0) <- next;
  Int.iter 1 (Array.length labels - 1) (fun i ->
    if labels.(i) = 0 then 0 else
    try IntMap.find nlabels labels.(i) with Not_found -> failwith "Labels.map: not found");
  labels

let remove_list labels l =
  let labels = Array.copy labels in
  Xlist.iter l (fun i -> labels.(i) <- 0);
  labels

let add_invisible_list labels l =
  let labels = Array.copy labels in
  Xlist.iter l (fun i ->   labels.(i) <- -1);
  labels
