def plot_tree(pdf_filename, clf):
    try:
        import pydot
    except ImportError, e:
        print("Please install 'pydot', it will require Graphviz.")
        sys.exit(1)

    from sklearn.externals.six import StringIO
    from sklearn.tree import export_graphviz

    dot_data = StringIO()
    export_graphviz(clf, out_file=dot_data)
    graph = pydot.graph_from_dot_data(dot_data.getvalue())
    graph.write_pdf(pdf_filename)

