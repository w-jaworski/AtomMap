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

open Big_int_Z

type 'a t

val empty : 'a t
val singleton : 'a list -> 'a list -> 'a t

val add : 'a t -> ('a list * 'a list) -> 'a t
val append : 'a t -> 'a t -> 'a t

val int_bucket_group : 'a t -> ('a -> int) -> ('a -> 'b) -> 'b t
val string_bucket_group : 'a t -> ('a -> string) -> ('a -> 'b) -> 'b t

val make_matchings : 'a t -> ('a * 'a) list list (* bylo wykomentowane *)

val fold : 'a t -> 'b -> ('b -> ('a list * 'a list) -> 'b) -> 'b

val list_fold : 'a list * 'a list -> 'b * 'b -> ('b -> 'a -> 'b) -> 'b * 'b
val list_map : 'a list * 'a list -> ('a -> 'b) -> 'b list * 'b list

val to_pair_list : 'a t -> ('a list * 'a list) list
val of_pair_list : ('a list * 'a list) list -> 'a t

val factorial : big_int -> int -> big_int

val count_matchings : 'a t -> big_int
(* val count_matchings : int t -> big_int *)
