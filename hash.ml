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

type hash_key_fun =  string array array ->
           Labels.t ->
           int ->
           Types.atom_props * (Types.bond * Types.atom_props) list -> string

(****
let key_of_valence = function
    Alifatic -> "X"
  | Aromatic -> "x"
****)
let key_of_bond = function
    Default -> "-"
  | Single -> "-"
  | Slash -> "-"
  | Backslash -> "-"
  | Double -> "="
  | Triple -> "#"
  | Aromatic -> "+"

let hash_key hash labels level (p,l) =
  if level = 0 then p.name ^ (try string_of_int (Labels.find labels p.id) with Not_found -> "") else
  let neigbours = String.concat "" (List.sort compare (Xlist.map l (fun (b,p) -> "-" ^ hash.(level-1).(p.id)))) in
  Printf.sprintf "{%sl%d%s}" p.name (Labels.get labels p.id) neigbours

let full_hash_key hash labels level (p,l) =
  if level = 0 then p.name ^ (try string_of_int (Labels.find labels p.id) with Not_found -> "") else
  let neigbours = String.concat "" (List.sort compare (Xlist.map l (fun (b,p) -> key_of_bond b ^ hash.(level-1).(p.id)))) in
  Printf.sprintf "{%s%dH%dF%dl%d%s}" p.name p.charge (Xlist.size p.hydrogens) (Xlist.size p.fluors)
    (Labels.get labels p.id) neigbours

let very_full_hash_key hash labels level (p,l) =
  if level = 0 then Printf.sprintf "{%s%dH%dF%dl%d}" p.name p.charge (Xlist.size p.hydrogens) (Xlist.size p.fluors) (Labels.get labels p.id) else
  let neigbours = String.concat "" (List.sort compare (Xlist.map l (fun (b,p) -> key_of_bond b ^ hash.(level-1).(p.id)))) in
  Printf.sprintf "{%s%dH%dF%dl%d%s}" p.name p.charge (Xlist.size p.hydrogens) (Xlist.size p.fluors)
    (Labels.get labels p.id) neigbours
(****
let unlabelled_hash_key hash graph labels level (p,l) =
  if level = 0 then Smiles.string_of_atom Smiles.NoId (Labels.make 1) p else
  let neigbours = String.concat "" (List.sort compare (Xlist.map l (fun (b,p) -> key_of_bond b ^ hash.(level-1).(p.id)))) in
  Printf.sprintf "{%s%s}" (Smiles.string_of_atom Smiles.NoId (Labels.make 1) p) neigbours
****)
let make_hash hash_key_fun level graph labels =
  let hash = Array.make (level+1) [| |] in
  Int.iter 0 level (fun i ->
    hash.(i) <- Array.map (fun (p,l) -> if p.id = 0 then "" else hash_key_fun hash labels i (p,l)) graph);
  hash
(****
let make_border_hash graph labels =
  Array.map (fun (p,l) ->
    if p.id = 0 then "" else
    if Labels.mem labels p.id then "" else
    if Xlist.fold l false (fun b (_,p) -> b || Labels.mem labels p.id) then
      let neigbours = String.concat "" (List.sort compare (Xlist.fold l [] (fun l (b,p) ->
        try ("-" ^ p.name ^ string_of_int (Labels.find labels p.id)) :: l
        with Not_found -> l))) in
      Printf.sprintf "%s%s" p.name neigbours
    else "") graph
****)
let make_border_hash2 level graph labels =
  let hash = Array.make (level+1) [| |] in
  Int.iter 0 (level-1) (fun i ->
    hash.(i) <- Array.map (fun (p,l) -> if p.id = 0 then "" else hash_key hash labels i (p,l)) graph);
  hash.(level) <- Array.map (fun (p,l) ->
    if p.id = 0 then "" else
    if Labels.mem labels p.id then "" else
    if Xlist.fold l false (fun b (_,p) -> b || Labels.mem labels p.id) then
      hash_key hash labels level (p,l)
    else "") graph;
  hash

let hash_bucket_group hash selection1 selection2 =
  Collection.string_pair_bucket_group selection1 selection2 (fun n -> if hash.(n) = "" then raise Not_found else hash.(n))

let hash_bucket_group2 hash1 hash2 selection1 selection2 =
(*   Collection.string_pair_bucket_group selection1 selection2 (fun n -> if hash1.(n)) (fun n -> if hash2.(n)) *)
  let map = Xlist.fold (Collection.to_list selection1) StringMap.empty (fun map n -> StringMap.add_inc map hash1.(n) ([n],[]) (fun (l,l2) -> n :: l,l2)) in
  Xlist.fold (Collection.to_list selection2) map (fun map n -> StringMap.add_inc map hash2.(n) ([],[n]) (fun (l,l2) -> l,n :: l2))

let hash_list_key hash selection =
  String.concat ";" (List.sort compare (Xlist.rev_map (Collection.to_list selection) (fun v -> hash.(v))))
