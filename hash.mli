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
open Types

type hash_key_fun =  hash array ->
           Labels.t ->
           int ->
           Types.atom_props * (Types.bond * Types.atom_props) list -> string

val key_of_bond : Types.bond -> string

val hash_key : hash_key_fun

val full_hash_key : hash_key_fun
val very_full_hash_key : hash_key_fun

val make_hash : hash_key_fun -> int -> graph -> Labels.t -> hash array

val make_border_hash2 : int -> graph -> Labels.t -> hash array

val hash_bucket_group : hash -> int Collection.collection -> int Collection.collection -> (int Collection.collection * int Collection.collection) Collection.collection
val hash_bucket_group2 : hash ->  hash -> int Collection.collection -> int Collection.collection -> (int list * int list) StringMap.t

val hash_list_key : hash -> int Collection.collection -> string
