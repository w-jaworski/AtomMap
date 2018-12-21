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
(****
type t

val get_connected_components : int list -> ('a * ('b * smile_atom_props) list) array -> IntSet.t -> int list list

(* val match_isomorphic_connected_components : int list list -> int list list -> (atom_props * ('a * atom_props) list) array -> IntSet.t -> IntSet.t -> string array array -> int -> int list Pair.t list *)

val isomorphic_partition : int -> 'a list list -> ('a * string array array * (smile_atom_props * ('b * smile_atom_props) list) array * IntSet.t) list -> 'a list list
(* val isomorphic_partition2 : int -> 'a list list -> ('a * (string array array * (atom_props * ('b * atom_props) list) array * IntSet.t) list) list -> 'a list list *)

val extend_matching3 : string array -> (smile_atom_props * ('a * smile_atom_props) list) array -> IntSet.t -> IntSet.t -> int list -> int list -> int Pair.t -> (int * int) list list
****)
exception Ambiguity of string
exception SolutionNotFound of string

val find_isomorphic_components : Types.graph ->
           int Collection.collection Collection.collection ->
           Labels.t ->
           int Collection.collection Collection.collection
           Collection.collection

(* etykietuje izomorficzne (zwn. hashe na zadanym poziomie) fragmenty cząsteczek uwzględniając dotychczasowe etykiety. *)
val label_reaction : int -> reaction -> Labels.t -> Labels.t list

(* następnie scala izomoriczne etykietowania poszczególnych cząsteczek*)
val disambiguate_labels : int -> reaction -> Labels.t list -> Labels.t list

val label_border_reaction : int -> reaction -> Labels.t -> Labels.t list

val match_reaction : int -> reaction -> Labels.t -> (int * int) Collection.collection Collection.collection (*(int list list * int list list) list*)
val match_reaction_full_key : int -> reaction -> Labels.t -> (int * int) Collection.collection Collection.collection (*(int list list * int list list) list*)

val disambiguate_matchings : reaction -> ((int * int) Collection.collection * (int * int) Collection.collection * IntSet.t) list ->
                                          ((int * int) Collection.collection * (int * int) Collection.collection * IntSet.t) list

(* val disambiguate_aroma : reaction list -> reaction list *)
val disambiguate_rearrangements : reaction list -> reaction list

type history
type action


val multilevel_label_reaction3 : Types.reaction ->
           history list ->
           Labels.t list -> action * history list * Labels.t list
