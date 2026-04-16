class dictWriter
// write out all values of a Dict verbatim as seen by Foam internally
// expands all #calc, #include etc. correctly. 
{
    protected: 
        word path_;
    public:
        word fileNameString_;
    // construct-and-write , custom dump name
    dictWriter(const IOdictionary& contentDict,const word fileNameString = "",  const word path = ".post/"):
        fileNameString_ (fileNameString),
        path_(path)
    {
        // if no name for file was given, default to dict name
        if  (fileNameString != "")  {fileNameString_ =fileNameString; Info << "in if" << endl;  }
        else {fileNameString_ = contentDict.name() ; }

        write(contentDict);
    }
    // constructor without write
    dictWriter(const word fileNameString, const word path = ".post/") :
        fileNameString_ (fileNameString),
        path_(path)
    {
    }

    const word& path(const word path = "")
    {
        if (path != "")
        {
            path_=path;
            mkDir (path_);
        }
        return path_;
    }
    // write call for a dict
    void write(const IOdictionary& contentDict)
    { 
        mkDir(path_);
        OFstream OFS (path_ + fileNameString_) ;
        writeHeader (contentDict,OFS); 
        // write without braces delimiters 
        forAllConstIter(dictionary, contentDict, iter)
        {
            OFS << iter() << nl;
        }

        OFS << endl; 
    }
 
    void writeHeader(const IOdictionary& contentDict, OFstream& OFS )
    {
        OFS  << "FoamFile\n"
            << "{\n"
            << "    version     2.0;\n"
            << "    format      ascii;\n"
            << "    class       dictionary;\n"
            << "    location    " << contentDict.filePath() << " ;\n"
            << "    object      " << contentDict.name() << ";\n"
            << "}\n"
            << "\n"; 
    }


};
 
