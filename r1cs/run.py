from prover import generate_proof
from verifier import verify_proof

solution_correct = 3
solution_incorrect = 4

proof_correct = generate_proof(3)
proof_incorrect = generate_proof(4)

print(verify_proof(proof_correct))
print(verify_proof(proof_incorrect))