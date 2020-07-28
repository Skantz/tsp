#load "str.cma";;

let query_node n node_lookup =
  (*let n : int = Random.int max_n in*)
  let x = (Hashtbl.find_opt node_lookup n) in
  match x with
  | None -> let () = Printf.printf "%d\n" n in
    let line = read_line () in
    let line = Str.split (Str.regexp_string " ") line |> List.map int_of_string in
    let () = Hashtbl.add node_lookup n (Some line) in
    line
  | Some x -> Option.get x
;;

let bfs start_node max_n node_lookup =
  let q = Queue.create () in
  let bound = int_of_float (1.0 /. Random.float 1.0) in
  let h = Hashtbl.create bound in
  let () = Queue.add start_node q in
  let i = ref 0 in
  while !i < bound do 
    let n = Queue.take q in
    let neighbors = query_node n node_lookup in
    List.iter (fun x -> Queue.add x q) neighbors;
    i := (!i + 1);
  done;
  (!i == bound)

let est_mst n eps max_w node_lookup =
  let res = Array.init max_w (fun x -> 0) in 
  let bound = int_of_float ((float_of_int max_w) /. (eps ** 2.0)) in
  let i = ref 0 in 
  while !i < bound do
    let start_node = Random.int n - 1 in
    let w = ref max_w in
    let bool_to_int t = if t then 1 else 0 in
    while !w > 0 do
      res.(!w - 1) <- res.(!w - 1) + bool_to_int (bfs start_node n node_lookup);
      w := !w + 1;
      ()
    done;
    ();
    i := !i + 1
  done;
  Array.map (fun x -> x * n / bound) res;

(*
let print_random_number max_n =
  let n = Random.int max_n in
  Printf.printf "%d\n" n
*)

Random.self_init ()

let n = read_int () 
let eps = read_float () -. 1.0 
let max_w = read_int ()

let node_lookup = Hashtbl.create 1000

print_random_number 15
let main () = query_node 10 node_lookup ()