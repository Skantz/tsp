#load "str.cma";;

let query_node n node_lookup =
  (*let n : int = Random.int max_n in*)
  let x = (Hashtbl.find_opt node_lookup n) in
  match x with
  | None -> let () = Printf.printf "%d\n" n in
    let line = read_line () in
    let line = Str.split (Str.regexp_string " ") line |> List.map int_of_string in
    let () = Hashtbl.add node_lookup n (Some line) in
    Some line
  | Some x -> x

let query_node n node_lookup = 

let bfs start_node max_n max_steps =
  let q = Queue.create () in
  (*let () = q.enqueue start_node in*)
  let bound = int_of_float (1.0 /. Random.float 1.0) in
  let h = Hashtbl.create bound in
  let i = 0 in
  let () = Queue.add start_node q in
  while i < bound do 
    let n = Queue.take q in
    let neighbors = query_node n in
    List.fold_left (fun x -> q.enqueue x) neighbors
    i := !i + 1
  done

let est_mst n eps max_w nbor_table =
  let bound = max_w /. (eps ** 2.0) in
  ()  

let print_random_number max_n =
  let n : int = Random.int max_n in
  Printf.printf "%d\n" n

Random.self_init()

let n = read_int ()
let eps = read_float () -. 1.0
let max_w = read_float ();;

let node_lookup = Hashtbl.create 1000;;


print_random_number 15;;
query_node 10 node_lookup;;