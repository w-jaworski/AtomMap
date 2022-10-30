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


type 'a collection

val empty : 'a collection

val of_list : 'a list -> 'a collection
val to_list : 'a collection -> 'a list

val size : 'a collection -> int

val sum : 'a collection -> 'a collection -> 'a collection

val partition : 'a collection -> ('a -> bool) -> 'a collection * 'a collection

val quotient : 'a collection -> ('a -> 'a -> bool) -> 'a collection collection

val reverse_projection : 'a collection collection -> 'a collection
val elem : 'a collection -> 'a
val singleton : 'a -> 'a collection
val is_empty : 'a collection -> bool
val flatten : 'a collection collection -> 'a collection

val annotate : 'a collection -> ('a -> 'b) -> ('a * 'b) collection

val deannotate : ('a * 'b) collection -> 'a collection

val combine : 'a collection -> 'b collection -> ('a -> 'b -> bool) -> ('a * 'b) collection

val string_pair_bucket_group : 'a collection -> 'a collection -> ('a -> string) -> ('a collection * 'a collection) collection
val string_pair_bucket_group2 : 'a collection -> 'a collection -> ('a -> string) -> ('a -> string) -> ('a collection * 'a collection) collection

val int_pair_bucket_group : 'a collection -> 'a collection -> ('a -> int) -> ('a collection * 'a collection) collection

val group : ('a * 'b) collection -> ('a * 'b collection) collection
val group2 : ('a * 'b) collection -> ('a collection * 'b) collection

val for_all : 'a collection -> ('a -> bool) -> bool
val exists : 'a collection -> ('a -> bool) -> bool

val flatten_map : 'a collection -> ('a -> 'b collection) -> 'b collection
val map : 'a collection -> ('a -> 'b) -> 'b collection

open Big_int_Z

val int_log : big_int -> int
val newton_symbol : int  -> int -> big_int

exception Empty

val number_of_matching_candidates : ('a collection * 'a collection) collection -> big_int
val number_of_exclusion_schemes : ('a collection * 'a collection) collection -> big_int

val generate_combinations : int -> 'a collection -> 'a collection collection
val generate_partition_combinations : int -> 'a collection -> ('a collection * 'a collection) collection

val generate_permutations : 'a collection * 'b collection -> ('a * 'b) collection collection

val number_of_product_permutations : ('a collection * 'b collection) collection -> big_int
val generate_product_permutations : ('a collection * 'b collection) collection -> ('a * 'b) collection collection

val number_of_products : 'a collection collection -> big_int
val product : 'a collection collection -> 'a collection collection
val flatten_product : 'a collection collection collection -> 'a collection collection

val get_exclusion_schemas_complete : ('a collection * 'a collection) collection -> 'a collection collection * 'a collection collection
val get_exclusion_schemas : ('a collection * 'a collection) collection -> 'a collection collection
val get_border_schemas : ('a collection * 'a collection) collection -> ('a collection * 'a collection) collection -> ('a collection * 'a collection) collection collection

val sort : 'a collection -> 'a collection

val uniq : 'a collection -> 'a collection

val intersection : 'a collection -> 'a collection -> 'a collection

val to_string : 'a collection -> string -> ('a -> string) -> string

val to_string_as_list : 'a collection -> ('a -> string) -> string

val print_cc : 'a collection collection -> string -> string -> ('a -> string) -> unit
val print_ccc : 'a collection collection collection -> string -> string -> string -> ('a -> string) -> unit
