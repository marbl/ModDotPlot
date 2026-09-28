from moddotplot.moddotplot import get_parser


def test_interactive_forward_defaults_to_false():
    args = get_parser().parse_args(["interactive", "--fasta", "sequence.fa"])

    assert args.forward is False


def test_interactive_delta_defaults_to_half_window():
    args = get_parser().parse_args(["interactive", "--fasta", "sequence.fa"])

    assert args.delta == 0.5


def test_interactive_accepts_forward_flag():
    args = get_parser().parse_args(
        ["interactive", "--fasta", "sequence.fa", "--forward"]
    )

    assert args.forward is True
