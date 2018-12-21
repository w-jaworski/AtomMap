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

let string_of_bond = function
    Default -> "-"
  | Single -> "-"
  | Slash -> "-"
  | Backslash -> "-"
  | Double -> "="
  | Triple -> "#"
  | Aromatic -> "+"

let string_of_charge = function
    0 -> ""
  | 1 -> "+"
  | -1 -> "-"
  | n -> if n > 0 then "+" ^ string_of_int n else "-" ^ string_of_int n

let make_hash graph i j =
  let p,l = graph.(i) in
  let charge = if p.charge = 0 then [] else [string_of_charge p.charge] in
  String.concat "" ([p.name] @ charge @ (List.sort compare (Xlist.fold l [] (fun found (bond,p) ->
    if p.id = j then found else
    (string_of_bond bond ^ p.name) :: found))))

let classify2 graph map (red,green) i =
   IntMap.fold map (red,green) (fun (red,green) j -> function
        1 -> (String.concat " " (List.sort compare [make_hash graph i j;make_hash graph j i])) :: red, green
      | 2 -> red, (String.concat " " (List.sort compare [make_hash graph i j;make_hash graph j i])) :: green
      | _ -> failwith "classify")

let classify (r,labels,center_bonds_modif,center_bonds,reid) =
  let redr,redp,greenr,greenp = IntMap.fold center_bonds ([],[],[],[]) (fun (redr,redp,greenr,greenp) i map ->
    if IntSet.mem r.reactant_ids i then
      let redr,greenr = classify2 r.graph map (redr,greenr) i in redr,redp,greenr,greenp
    else let redp,greenp = classify2 r.graph map (redp,greenp) i in redr,redp,greenr,greenp) in
  redr,redp,greenr,greenp

let count_red_oh (r,labels,center_bonds_modif,center_bonds,reid) =
  IntMap.fold center_bonds 0 (fun no_red_oh i map ->
    if IntSet.mem r.reactant_ids i then
      IntMap.fold map no_red_oh (fun no_red_oh j -> function
        1 ->
          let i,j = if (fst r.graph.(i)).name = "O" then j,i else i,j in
          if (fst r.graph.(i)).name = "C" && make_hash r.graph j i = "O-H" then no_red_oh + 1 else no_red_oh
      | _ -> no_red_oh)
    else no_red_oh)

let to_string (redr,redp,greenr,greenp) =
  "red reactants " ^ String.concat " " (List.sort compare redr) ^ " " ^
  "green reactants " ^ String.concat " " (List.sort compare greenr) ^ " " ^
  "red products " ^ String.concat " " (List.sort compare redp) ^ " " ^
  "green products " ^ String.concat " " (List.sort compare greenp)

let disambiguate solutions =
  let solutions = Collection.of_list solutions in
  let solutions_ann = Collection.map solutions (fun solution -> to_string (classify solution), solution) in
  let grouped_solution = Collection.map (Collection.group solutions_ann) (fun (_,sols) -> sols) in
  Collection.to_list (Collection.reverse_projection grouped_solution)

let selection_schemata = [
  ["C-C=O O-H"; "H O-C"];
  ["N=O O-H"; "H O-C"];
  ["C-C=C O-H"; "H O-C"];
  ]

let rec check_schema redr = function
    [] -> true
  | a :: l -> if Xlist.mem redr a then check_schema redr l else false

let select solutions =
  if Xlist.size solutions < 2 then solutions else
  let selected,rest = Xlist.fold solutions ([],[]) (fun (selected,rest) solution ->
(*     printf "solution %d\n" (count_red_oh solution); *)
    let redr,_,_,_ = classify solution in
(*     Xlist.iter redr (printf "\"%s\"\n"); *)
    let b = Xlist.fold selection_schemata false (fun b l -> b || check_schema redr l) in
    if b then solution :: selected, rest else selected, solution :: rest) in
  let solutions = if selected = [] then rest else selected in
  fst (Xlist.fold solutions ([],max_int) (fun (selected,n) solution ->
    let m = count_red_oh solution in
    if m < n then [solution],m else
    if m > n then selected,n else
    solution :: selected, m))

let rec is_cycle_rec n g i visited found =
  if n = 0 then found else
  let visited = if n = 1 then IntSet.remove visited i else visited in
  let visited = IntSet.union visited found in
  let found = IntSet.fold found IntSet.empty (fun found i ->
    let set = try IntMap.find g i with Not_found -> IntSet.empty in
    IntSet.fold set found (fun found i ->
      if IntSet.mem visited i then found else IntSet.add found i)) in
  is_cycle_rec (n-1) g i visited found

let is_cycle g i j =
  let set = is_cycle_rec 5 g i (IntSet.singleton i) (IntSet.singleton j) in
  IntSet.mem set i

let remove_cycles solutions =
  Xlist.map solutions (fun (r,labels,center_bonds_modif,center_bonds,reid) ->
    let g = IntMap.fold center_bonds IntMap.empty (fun g i map ->
      IntMap.fold map g (fun g j v ->
        if v = 2 then
          let g = IntMap.add_inc g i (IntSet.singleton j) (fun set -> IntSet.add set j) in
          IntMap.add_inc g j (IntSet.singleton i) (fun set -> IntSet.add set i)
        else g)) in
    let center_bonds = IntMap.fold center_bonds IntMap.empty (fun center_bonds i map ->
      let map = IntMap.fold map IntMap.empty (fun map j v ->
        if v = 2 then if is_cycle g i j then map else IntMap.add map j v
        else IntMap.add map j v) in
      if IntMap.is_empty map then center_bonds else IntMap.add center_bonds i map) in
    r,labels,center_bonds_modif,center_bonds,reid)
